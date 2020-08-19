import shutil
from os import path
from pathlib import Path

from embers.sat_utils.sat_ephemeris import (ephem_data, epoch_ranges,
                                            epoch_time_array, load_tle,
                                            sat_pass, sat_plot, save_ephem)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


def test_load_tle_sats():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    sat_id = sats[0].model.satnum
    assert sat_id == 25986


def test_load_tle_epochs():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    assert epochs[0] == 2458738.5


def test_epoch_ranges_length():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    assert len(epoch_range) == 311


def test_epoch_time_array_index():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    assert index_epoch == 0


def test_epoch_time_array_arr():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    assert type(t_arr).__name__ == "Time"


def test_sat_pass_passes_1():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    assert data[0][0][0] == 417


def test_sat_pass_passes_2():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=1, cadence=10)
    data = sat_pass(sats, t_arr, 1, location=(-26.703319, 116.670815, 337.83))
    assert data[0][0][0] == 0


def test_sat_pass_alt_az():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    assert type(data[1]) == type(data[2])


def test_sat_pass_alt_err():
    tle_file = f"{test_data}/sat_utils/TLE/44387.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    assert data is None


def test_ephem_data_time():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    passes, alt, az = data
    time_array, sat_alt, sat_az = ephem_data(t_arr, passes[0], alt, az)
    assert time_array.shape[0] == 98


def test_ephem_data_alt():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    passes, alt, az = data
    time_array, sat_alt, sat_az = ephem_data(t_arr, passes[0], alt, az)
    assert sat_alt.shape[0] == 98


def test_ephem_data_az():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    passes, alt, az = data
    time_array, sat_alt, sat_az = ephem_data(t_arr, passes[0], alt, az)
    assert sat_az.shape[0] == 98


def test_sat_plot():
    tle_file = f"{test_data}/sat_utils/TLE/25986.txt"
    sats, epochs = load_tle(tle_file)
    epoch_range = epoch_ranges(epochs)
    t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch=0, cadence=10)
    data = sat_pass(sats, t_arr, 0, location=(-26.703319, 116.670815, 337.83))
    passes, alt, az = data
    time_array, sat_alt, sat_az = ephem_data(t_arr, passes[0], alt, az)
    plt = sat_plot(25986, sat_alt, sat_az)
    assert plt.__name__ == "matplotlib.pyplot"


def test_save_ephem_empty_file():
    error = save_ephem(
        12345,
        f"{test_data}/sat_utils/TLE",
        10,
        (-26.703319, 116.670815, 337.83),
        0.5,
        f"{test_data}/sat_utils/ephem_tmp",
    )
    assert error == f"File {test_data}/sat_utils/TLE/12345 is empty, skipping"


def test_save_ephem_plot():
    save_ephem(
        44387,
        f"{test_data}/sat_utils/TLE",
        10,
        (-26.703319, 116.670815, 337.83),
        0.5,
        f"{test_data}/sat_utils/ephem_tmp",
    )
    png = Path(f"{test_data}/sat_utils/ephem_tmp/ephem_plots/44387.png")
    assert png.is_file() is True
    if png.is_file() is True:
        shutil.rmtree(f"{test_data}/sat_utils/ephem_tmp")
