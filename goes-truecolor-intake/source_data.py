#!/usr/bin/env python
# coding: utf-8


import datetime
from pathlib import Path

import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import numpy as np
import rioxarray as rxr
import satdata
import satpy
import xarray as xr
from loguru import logger
from suntime import Sun
from tqdm.auto import tqdm

DEBUG = True
DATETIME_FORMAT = "%Y%m%dT%H%M%SZ"


def plot_domain(bbox):
    fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))

    ax.set_extent(bbox, crs=ccrs.PlateCarree())
    ax.gridlines(draw_labels=["left", "bottom"])
    ax.coastlines()


def calc_minmax_time_for_bbox(bbox, date):
    west_edge = [(bbox[2] + bbox[3]) / 2, bbox[0]]
    east_edge = [(bbox[2] + bbox[3]) / 2, bbox[1]]

    t_min = Sun(*east_edge).get_sunrise_time(date)
    t_max = Sun(*west_edge).get_sunset_time(date)

    return t_min, t_max


def _parse_time_from_source_filename(filename):
    return satdata.Goes16AWS.parse_key(filename, parse_times=True)["start_time"]


def download_radiance_source_files(
    bbox, t_min, t_max, data_path_root=Path("source_data/")
):
    """ """
    select_nearest_to_t_min = False
    if t_min == t_max:
        dt_total = datetime.timedelta(minutes=30)
        select_nearest_to_t_min = True
    else:
        dt_total = t_max - t_min
    t_center = t_min + dt_total / 2

    cli = satdata.Goes16AWS(local_storage_dir=data_path_root)
    files_per_channel = {}
    for channel in [1, 2, 3]:
        keys = cli.query(
            # XXX
            time=t_center.replace(tzinfo=None),
            region="F",
            debug=True,
            dt_max=dt_total / 2,
            channel=channel,
        )
        if DEBUG:
            N_keys = len(keys) // 2
            keys = keys[N_keys - 1 : N_keys + 2]
        filenames = cli.download(
            keys,
        )
        files_per_channel[channel] = filenames

    if select_nearest_to_t_min:
        scene_times = np.array(
            [
                _parse_time_from_source_filename(filename=fn)
                for fn in files_per_channel[0]
            ]
        )
        t_dist = np.abs(scene_times - np.datetime64(t_min))
        idx_nearest = np.argmin(t_dist)

        files_per_channel = {
            c: fns[idx_nearest] for (c, fns) in files_per_channel.items()
        }

    return files_per_channel


def _merge_multichannel_sources(files_per_channel):
    channel_files_by_timestamp = {}
    N_channels = len(files_per_channel)
    for channel, channel_files in files_per_channel.items():
        for ch_filename in channel_files:
            file_timestamp = _parse_time_from_source_filename(filename=ch_filename)
            time_group = channel_files_by_timestamp.setdefault(file_timestamp, {})
            time_group[channel] = ch_filename

    scene_filesets = []

    for timestamp in sorted(channel_files_by_timestamp.keys()):
        timestamp_files = channel_files_by_timestamp[timestamp]

        if len(timestamp_files) == N_channels:
            scene_filesets.append(
                [timestamp_files[channel] for channel in files_per_channel.keys()]
            )
        else:
            logger.warn(
                f"Only {len(timestamp_files)} were found for timestamp {timestamp}"
                " so this timestamp will be excluded"
            )

    return scene_filesets


def create_truecolor_scene_files(t_min, t_max, bbox, make_netcdf=False, data_path="."):
    files_per_channel = download_radiance_source_files(
        t_min=t_min,
        t_max=t_max,
        bbox=bbox,
        data_path_root=Path(data_path) / "source_data",
    )

    scene_filesets = _merge_multichannel_sources(files_per_channel=files_per_channel)

    for scene_files in tqdm(scene_filesets):
        start_time = _parse_time_from_source_filename(scene_files[0])
        start_time_str = start_time.strftime(DATETIME_FORMAT)
        geotiff_filepath = Path(f"{start_time_str}.tif")
        if not geotiff_filepath.exists():
            create_cropped_truecolor_geotiff(
                scene_filepaths=scene_files,
                geotiff_filepath=geotiff_filepath,
                bbox=bbox,
            )

        if make_netcdf:
            geotiff_to_netcdf(
                geotiff_filepath=geotiff_filepath, channel1_filepath=scene_files[0]
            )


def create_cropped_truecolor_geotiff(scene_filepaths, geotiff_filepath, bbox):
    """
    Load radiance files (assumed to be channels 1, 2 and 3) in
    `scene_filepaths` and cropped with bounding-box (`bbox` given as list in
    WESN format) create a true-color RGB composite geotiff and save to
    `geotiff_filepath`
    """
    # Creating truecolor composites
    scene = satpy.Scene(scene_filepaths, reader="abi_l1b")

    # composite
    scene.load(["true_color"])

    # satpy bbox: (xmin, ymin, xmax, ymax), WSEN
    ll_bbox = [bbox[0], bbox[2], bbox[1], bbox[3]]
    scene_cropped = scene.crop(ll_bbox=ll_bbox)

    # it is necessary to "resample" here because the different channels are at
    # different spatial resolution. By not passing in an "area" the highest
    # resolution possible will be used
    new_scn = scene_cropped.resample(resampler="native")

    new_scn.save_dataset("true_color", writer="geotiff", filename=str(geotiff_filepath))


def geotiff_to_netcdf(geotiff_filepath, channel1_filepath):
    """
    Convert geotiff image in filepath `geotiff_filepath` to a netCDF file that
    can be plotted with `matplotlib.imshow` which contains projection
    information in a CF-compliant manner (so that it can be read by `metpy`).

    NB: Because of the incompatibility between cartopy and pyproj at the moment we
    can't read the projection info directly from the geotiff, but instead must
    get it from one of the radiance channel source files. This projection
    information will be read from `channel1_filepath`

    """
    # Convert geotiff into netCDF file with RGB image data
    ds_ch1 = xr.open_dataset(channel1_filepath)

    # read tiff with rioxarray
    da = rxr.open_rasterio(geotiff_filepath)

    # remove attribute that we can't easily serialise and scale from uint8 in
    # range 0..255 to float in range 0...1
    da_img_floats = da.copy().drop("spatial_ref").expand_dims("t") / 255

    # copy over projection info in a CF-compliant manner so that metpy can parse it
    grid_mapping = ds_ch1.Rad.grid_mapping
    da_img_floats.coords[grid_mapping] = ds_ch1[grid_mapping].drop(
        ["t", "x_image", "y_image"]
    )
    da_img_floats.attrs["grid_mapping"] = grid_mapping
    da_img_floats.name = "true_color"

    netcdf_filepath = geotiff_filepath.parent / geotiff_filepath.name.replace(
        ".tif", ".nc"
    )
    da_img_floats.to_netcdf(netcdf_filepath)
