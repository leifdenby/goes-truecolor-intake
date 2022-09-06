#!/usr/bin/env python
# coding: utf-8
import isodate

from .source_data import calc_minmax_time_for_bbox, create_truecolor_scene_files

def _bbox_arg(s):
    bbox = [float(v) for v in s.split(',')]

    if len(bbox) != 4:
        raise Exception("`bbox` should be a list of [lat_S, lon_W, lat_N, lon_E]")
    return bbox


def main():
    import argparse

    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "date",
        type=isodate.parse_date,
        help="Date to download data for, for example 2020-01-01",
    )
    argparser.add_argument(
        "bbox",
        type=_bbox_arg,
        help="Domain bounding-box W,E,S,N for example -60,-50,11,16",
    )
    args = argparser.parse_args()

    bbox = args.bbox
    date = args.date
    t_min, t_max = calc_minmax_time_for_bbox(date=date, bbox=bbox)

    create_truecolor_scene_files(t_min=t_min, t_max=t_max, bbox=bbox)


if __name__ == "__main__":
    main()
