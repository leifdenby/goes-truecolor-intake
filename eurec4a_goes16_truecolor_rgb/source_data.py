#!/usr/bin/env python
# coding: utf-8


import datetime
import satdata
import satpy
from pathlib import Path
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import rioxarray as rxr
from satpy.writers import get_enhanced_image


from suntime import Sun, SunTimeException
import datetime


DOMAIN_BBOX = [-60, -50, 11, 16]  # WESN


def plot_domain():
    fig, ax = plt.subplots(subplot_kw=dict(projection=ccrs.PlateCarree()))

    ax.set_extent(DOMAIN_BBOX, crs=ccrs.PlateCarree())
    ax.gridlines(draw_labels=["left", "bottom"])
    ax.coastlines()


def calc_minmax_time_for_bbox(bbox, date):
    west_edge = [(bbox[2] + bbox[3]) / 2, bbox[0]]
    east_edge = [(bbox[2] + bbox[3]) / 2, bbox[1]]

    t_min = Sun(*east_edge).get_sunrise_time(date)
    t_max = Sun(*west_edge).get_sunset_time(date)

    return t_min, t_max


def download_daytime_radiance_source_files(
    date, bbox, data_path_root=Path("source_data/")
):
    """ """

    t_min, t_max = calc_minmax_time_for_bbox(date=date, bbox=bbox)
    dt_daytime = t_max - t_min

    t_center = t_min + dt_daytime / 2

    cli = satdata.Goes16AWS()
    filenames = []
    for channel in [1, 2, 3]:
        keys = cli.query(
            # XXX
            time=t_center.replace(tzinfo=None),
            region="F",
            debug=True,
            dt_max=dt_daytime / 2,
            channel=channel,
        )
        fn = cli.download(keys)
        filenames.append(fn)

    return filenames


def main():
    date = datetime.datetime.now()
    print(date)
    bbox = DOMAIN_BBOX
    files = download_daytime_radiance_source_files(date=date, bbox=bbox)
    print(files)
    # scene_meta = satdata.Goes16AWS.parse_key(scene_filepaths[0], parse_times=True)
    # start_time = scene_meta["start_time"]

    # # make ISO8601 identifier
    # start_time_str = start_time.isoformat().replace("-", "").replace(":", "") + "Z"
    # filename = f"GOES16_true_color_{start_time_str}.tif"


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
    ll_bbox = [DOMAIN_BBOX[0], DOMAIN_BBOX[2], DOMAIN_BBOX[1], DOMAIN_BBOX[3]]

    scene_cropped = scene.crop(ll_bbox=ll_bbox)

    # it is necessary to "resample" here because the different channels are at
    # different spatial resolution. By not passing in an "area" the highest
    # resolution possible will be used
    new_scn = scene_cropped.resample(resampler="native")

    new_scn.save_dataset("true_color", writer="geotiff", filename=geotiff_filepath)


def geotiff(geotiff_filepath):
    """ """
    # Convert geotiff into netCDF file with RGB image data
    ds_ch1 = xr.open_dataset(filenames[0])

    # read tiff with rioxarray
    da = rxr.open_rasterio(filename)

    # remove attribute that we can't easily serialise
    da2 = da.copy().drop("spatial_ref").expand_dims("t") / 255

    # copy over projection info in a CF-compliant manner so that metpy can parse it
    grid_mapping = ds_ch1.Rad.grid_mapping
    da2.coords[grid_mapping] = ds_ch1[grid_mapping].drop(["t", "x_image", "y_image"])
    da2.attrs["grid_mapping"] = grid_mapping
    da2.name = "true_color"

    da2.to_netcdf(filename.replace(".tif", ".nc"))


if __name__ == "__main__":
    main()
