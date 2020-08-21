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


def test_download_meta():
    download_meta(
        "2019-10-01", "2019-10-10", 1, f"{test_data}/mwa_utils/mwa_meta_tmp", 0
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


def test_point_integration_point_0():
    pointings, int_hours = point_integration(f"{test_data}/mwa_utils")
    assert pointings[0] == 0


def test_point_integration_int_0():
    pointings, int_hours = point_integration(f"{test_data}/mwa_utils")
    assert round(int_hours[0]) == 24


def test_pointing_hist():
    pointings, int_hours = point_integration(f"{test_data}/mwa_utils")
    pointing_hist(pointings, int_hours, 4, f"{test_data}/mwa_utils/tmp")
    hist = Path(f"{test_data}/mwa_utils/tmp/pointing_integration.png")
    assert hist.is_file()
    if hist.is_file():
        shutil.rmtree(f"{test_data}/mwa_utils/tmp")


def test_rf_obs_times_start():
    obs_time, obs_gps, obs_gps_end = rf_obs_times(
        "2019-10-01", "2019-10-01", "Australia/Perth"
    )
    assert obs_time[0] == "2019-10-01-00:00"


def test_rf_obs_times_gps():
    obs_time, obs_gps, obs_gps_end = rf_obs_times(
        "2019-10-01", "2019-10-01", "Australia/Perth"
    )
    assert obs_gps[0] == 1253894418.0


def test_rf_obs_times_interval():
    obs_time, obs_gps, obs_gps_end = rf_obs_times(
        "2019-10-01", "2019-10-01", "Australia/Perth"
    )
    assert obs_gps_end[0] - obs_gps[0] == 1800


def test_obs_pointings():
    obs_pointings(
        "2019-10-01", "2019-10-01", "Australia/Perth", f"{test_data}/mwa_utils"
    )
    obs_p = Path(f"{test_data}/mwa_utils/obs_pointings.json")
    assert obs_p.is_file()
    if obs_p.is_file():
        obs_p.unlink()


def test_tile_integration():
    obs_pointings(
        "2019-10-01", "2019-10-01", "Australia/Perth", f"{test_data}/mwa_utils"
    )
    tile_ints = tile_integration(
        f"{test_data}/mwa_utils", f"{test_data}/rf_tools/rf_data"
    )
    assert tile_ints["rf0XX"] == [0, 0, 0, 0]
    obs_p = Path(f"{test_data}/mwa_utils/obs_pointings.json")
    if obs_p.is_file():
        obs_p.unlink()


def test_plt_hist_array():
    obs_pointings(
        "2019-10-01", "2019-10-01", "Australia/Perth", f"{test_data}/mwa_utils"
    )
    tile_ints = tile_integration(
        f"{test_data}/mwa_utils", f"{test_data}/rf_tools/rf_data"
    )
    plt_hist_array(tile_ints, f"{test_data}/mwa_utils")
    obs_p = Path(f"{test_data}/mwa_utils/obs_pointings.json")
    if obs_p.is_file():
        obs_p.unlink()
    tile_hist = Path(f"{test_data}/mwa_utils/tiles_pointing_integration.png")
    assert tile_hist.is_file()
    if tile_hist.is_file():
        tile_hist.unlink()


def test_mwa_point_meta():
    mwa_point_meta(
        "2019-10-01",
        "2019-10-10",
        1,
        4,
        "Australia/Perth",
        f"{test_data}/rf_tools/rf_data",
        f"{test_data}/mwa_utils/mwa_meta_tmp",
        wait=0,
    )
    assert len(list(Path(f"{test_data}/mwa_utils/mwa_meta_tmp").glob("*"))) == 5
    shutil.rmtree(f"{test_data}/mwa_utils/mwa_meta_tmp")
