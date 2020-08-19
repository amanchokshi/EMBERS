from os import path
from pathlib import Path

from embers.sat_utils.sat_ephemeris import (ephem_data, epoch_ranges,
                                            epoch_time_array, load_tle,
                                            sat_pass, sat_plot, save_ephem)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


tle_file = f"{test_data}/TLE/25986.txt"

sats, epochs = load_tle(tle_file)
print(sats[0].model.satnum)
_, epochs = load_tle(tle_file)
epoch_range = epoch_ranges(epochs)
print(len(epoch_range))


def test_load_tle_sats():
    sats, _ = load_tle(tle_file)
    sat_id = sats[0].model.satnum
    assert sat_id == 25986


def test_load_tle_epochs():
    _, epochs = load_tle(tle_file)
    assert epochs[0] == 2458756.5


def test_epoch_ranges_length():
    _, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    assert len(epoch_range) == 20
