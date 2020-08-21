import shutil
from os import path
from pathlib import Path

from embers.mwa_utils.mwa_pointings import (clean_meta_json, combine_pointings,
                                            download_meta, mwa_point_meta,
                                            obs_pointings, plt_hist_array,
                                            point_integration, pointing_hist,
                                            rf_obs_times, tile_integration)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))

start_gps, stop_gps, obs_length, pointings = clean_meta_json(f"{test_data}/mwa_utils/")


def test_download_meta():
    download_meta(
        "2019-10-01", "2019-10-10", 1, f"{test_data}/mwa_utils/mwa_meta_tmp", wait=0
    )
    meta = Path(f"{test_data}/mwa_utils/mwa_meta_tmp/mwa_pointings/page_001.json")
    assert meta.is_file()
    if meta.is_file():
        shutil.rmtree(f"{test_data}/mwa_utils/mwa_meta_tmp")


def test_clean_meta_json_start():
    start_gps, stop_gps, obs_length, pointings = clean_meta_json(
        f"{test_data}/mwa_utils"
    )
    assert start_gps[0] == 1253960760


def test_clean_meta_json_point():
    start_gps, stop_gps, obs_length, pointings = clean_meta_json(
        f"{test_data}/mwa_utils"
    )
    assert pointings[0] == 0


def test_combine_pointings():
    start_gps, stop_gps, obs_length, pointings = clean_meta_json(
        f"{test_data}/mwa_utils"
    )
    combine_pointings(
        start_gps, stop_gps, obs_length, pointings, f"{test_data}/mwa_utils/tmp"
    )
    mwa_meta = Path(f"{test_data}/mwa_utils//tmp/mwa_pointings.json")
    assert mwa_meta.is_file()
    if mwa_meta.is_file():
        mwa_meta.unlink()
