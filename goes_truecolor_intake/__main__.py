#!/usr/bin/env python
# coding: utf-8
import isodate

from .s3_upload import create_catalog, upload_files
from .source_data import calc_minmax_time_for_bbox, create_truecolor_scene_files


def _bbox_arg(s):
    bbox = [float(v) for v in s.split(",")]

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
    argparser.add_argument("--local-path", default=".")
    argparser.add_argument("--s3-access-key", default=None)
    argparser.add_argument("--s3-secret-key", default=None)
    argparser.add_argument("--s3-bucket-name", default=None)
    argparser.add_argument("--target-hostname", default="local")
    args = argparser.parse_args()

    if args.target_hostname != "local":
        reqd_args = "s3_access_key s3_secret_key s3_bucket_name".split()
        if any(getattr(args, v) is None for v in reqd_args):
            raise Exception(
                "To upload to an S3 host you need to provide: " + ", ".join(reqd_args)
            )

    local_path = args.local_path
    name = "goesrgb"

    bbox = args.bbox
    date = args.date
    t_min, t_max = calc_minmax_time_for_bbox(date=date, bbox=bbox)

    create_truecolor_scene_files(
        t_min=t_min, t_max=t_max, bbox=bbox, data_path=local_path
    )

    create_catalog(
        catalog_path=local_path,
        name=name,
        host_url=args.target_hostname,
        remote_path=args.s3_bucket_name,
    )

    if args.target_hostname != "local":
        upload_files(
            local_path=local_path,
            host_url=args.target_hostname,
            remote_path=args.s3_bucket_name,
            s3_access_key=args.s3_access_key,
            s3_secret_key=args.s3_secret_key,
        )


if __name__ == "__main__":
    main()
