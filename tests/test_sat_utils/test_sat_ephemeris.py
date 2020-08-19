from os import path

from embers.sat_utils.sat_ephemeris import (ephem_data, epoch_ranges,
                                            epoch_time_array, load_tle,
                                            sat_pass, sat_plot, save_ephem)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


tle_file = f"{test_data}/TLE/25986.txt"

sats, epochs = load_tle(tle_file)
print(epochs)
epoch_range = epoch_ranges(epochs)
t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
print(data)


def test_load_tle_sats():
    sat_id = sats[0].model.satnum
    assert sat_id == 25986


def test_load_tle_epochs():
    assert epochs[0] == 2458738.5


def test_epoch_ranges_length():
    assert len(epoch_range) == 311


def test_epoch_time_array_index():
    assert index_epoch == 0


def test_epoch_time_array_arr():
    assert type(t_arr).__name__ == "Time"


def test_sat_pass_passes_1():
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    assert data[0][0][0] == 417


def test_sat_pass_passes_2():
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=1, cadence=10)
    data = sat_pass(sats, t_arr, 1, location=(-26.703319, 116.670815, 337.83))
    assert data[0][0][0] == 0


def test_sat_pass_alt_az():
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    assert type(data[1]) == type(data[2])


def test_sat_pass_alt_err():
    tle_file = f"{test_data}/TLE/44387.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    assert data is None
