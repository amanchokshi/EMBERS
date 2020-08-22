from os import path
from pathlib import Path

from embers.tile_maps.tile_maps import (check_pointing, mwa_clean_maps,
                                        plt_channel, plt_clean_maps,
                                        plt_fee_fit, plt_sat_maps,
                                        project_tile_healpix,
                                        rf_apply_thresholds, rfe_batch_cali,
                                        rfe_calibration, rfe_collate_cali,
                                        tile_maps_batch)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))

nside = 32


def test_check_pointing_0():
    point = check_pointing(
        "2019-10-01-23:30", f"{test_data}/tile_maps/obs_pointings.json"
    )
    assert point == 0


def test_check_pointing_2():
    point = check_pointing(
        "2019-10-01-23:00", f"{test_data}/tile_maps/obs_pointings.json"
    )
    assert point == 2


def test_check_pointing_4():
    point = check_pointing(
        "2019-10-02-00:00", f"{test_data}/tile_maps/obs_pointings.json"
    )
    assert point == 4


def test_check_pointing_41():
    point = check_pointing(
        "2019-10-07-05:30", f"{test_data}/tile_maps/obs_pointings.json"
    )
    assert point == 41


def test_check_pointing_none():
    point = check_pointing(
        "2019-10-08-00:00", f"{test_data}/tile_maps/obs_pointings.json"
    )
    assert point is None
