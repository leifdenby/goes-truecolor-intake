#!/usr/bin/env python
# coding: utf-8

from pathlib import Path

import requests
import s3fs
import yaml
from intake import open_catalog
from loguru import logger

from .source_data import DATETIME_FORMAT

requests.packages.urllib3.util.ssl_.DEFAULT_CIPHERS = "AES256-SHA"  # noqa

LOCAL_CATALOG_FILENAME = "catalog.yml"

BACKEND_SLEEP_TIME = 10


def upload_files(local_path, host_url, remote_path, s3_access_key, s3_secret_key):
    logger.info(f"storing in {host_url} {remote_path}...", end="", flush=True)

    endpoint_url = f"https://{host_url}"
    logger.info(f"Uploading files to remote S3 bucket {endpoint_url}")

    fs = s3fs.S3FileSystem(
        anon=False,
        key=s3_access_key,
        secret=s3_secret_key,
        client_kwargs={"endpoint_url": endpoint_url},
    )

    files_to_copy = list(Path(local_path).glob("*.tif"))

    logger.info(fs.ls("eurec4a-environment"))

    for fp_local in files_to_copy:
        fs.put(str(fp_local), str(Path("eurec4a-environment/goesrgb") / fp_local.name))


def create_catalog(catalog_path, name, host_url, remote_path):
    cat = {}
    cat["sources"] = {}

    if host_url == "local":
        url_root = "{{CATALOG_DIR}}"
    else:
        url_root = f"s3://{remote_path}/{name}"

    cat_entry = {
        "description": "Truecolor RGB of GOES-16",
        "driver": "rasterio",
        "args": {
            "urlpath": url_root + "/{time:" + DATETIME_FORMAT + "}.tif",
            "engine": "rioxarray",
            "concat_dim": "time",
            "chunks": {"band": 1, "x": 100, "y": 100},
        },
    }

    if host_url != "local":
        endpoint_url = f"https://{host_url}"
        cat_entry["args"]["storage_options"] = {
            "anon": True,
            "client_kwargs": {"endpoint_url": endpoint_url},
        }
        cat_entry["cache"] = [
            {
                "argkey": "urlpath",
                "regex": "goesrgb",
                "type": "file",
            }
        ]

    cat["sources"][name] = cat_entry

    with open(Path(catalog_path) / LOCAL_CATALOG_FILENAME, "w") as fh:
        fh.write(yaml.dump(cat, default_flow_style=None))

    logger.info(f"Wrote intake catalog file to `{LOCAL_CATALOG_FILENAME}`")


def _test_local_catalog_entry(name):
    cat = open_catalog(LOCAL_CATALOG_FILENAME)
    print(cat[name].to_dask())
