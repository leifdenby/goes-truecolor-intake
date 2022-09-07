#!/usr/bin/env python
# coding: utf-8

import requests
import s3fs
import yaml
from intake import open_catalog
from loguru import logger

from .source_data import DATETIME_FORMAT

requests.packages.urllib3.util.ssl_.DEFAULT_CIPHERS = "AES256-SHA"  # noqa

LOCAL_CATALOG_FILENAME = "catalog.yml"

BACKEND_SLEEP_TIME = 10


def _upload(ds, host_url, remote_path, auth_info):
    logger.info(f"storing in {host_url} {remote_path}...", end="", flush=True)
    s3_access_key = auth_info["s3_access_key"]
    s3_secret_key = auth_info["s3_secret_key"]

    fs = s3fs.S3FileSystem(
        anon=False,
        key=s3_access_key,
        secret=s3_secret_key,
        client_kwargs={"endpoint_url": host_url},
    )

    mapper = fs.get_mapper(remote_path)
    ds.to_zarr(mapper, consolidated=True)


def create_catalog(name, host_url, remote_path):
    cat = {}
    cat["sources"] = {}

    cat_entry = {
        "description": "Truecolor RGB of GOES-16",
        "driver": "rasterio",
        "args": {
            # "urlpath": f"simplecache::s3://{remote_path}/{name}/\{{DATETIME_FORMAT:time}\}.tif",
            "urlpath": "{{CATALOG_DIR}}/{time:" + DATETIME_FORMAT + "}.tif",
            "storage_options": {
                "s3": {
                    "anon": True,
                    "client_kwargs": {"endpoint_url": host_url},  # noqa
                }
            },
            "engine": "rioxarray",
            "concat_dim": "time",
            "chunks": {"band": 1, "x": 100, "y": 100},
        },
    }
    cat["sources"][name] = cat_entry

    with open(LOCAL_CATALOG_FILENAME, "w") as fh:
        fh.write(yaml.dump(cat, default_flow_style=None))

    logger.info(f"Wrote intake catalog file to `{LOCAL_CATALOG_FILENAME}`")


def _test_local_catalog_entry(name):
    cat = open_catalog(LOCAL_CATALOG_FILENAME)
    print(cat[name].to_dask())
