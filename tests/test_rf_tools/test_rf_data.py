import pytest
from embers.rf_tools.rf_data import (
    read_data,
    tile_names,
    tile_pairs,
    time_tree,
    plt_waterfall,
    single_waterfall,
    batch_waterfall,
)


def test_read_data_power_shape():
    power, _ = read_data(
        rf_file="../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert power.shape == (16655, 112)


def test_read_data_times_shape():
    _, times = read_data(
        rf_file="../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert list(times.shape)[0] == 16655


def test_read_data_power_times_shape():
    power, times = read_data(
        rf_file="../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert list(times.shape)[0] == list(power.shape)[0]


def test_read_data_return_type():
    _, times = read_data(
        rf_file="../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert type(times).__name__ == "ndarray"


def test_read_data_time_arrow():
    _, times = read_data(
        rf_file="../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert times[0] <= times[-1]


def test_tile_names_first():
    tiles = tile_names()
    assert tiles[0] == "rf0XX"


def test_tile_names_last():
    tiles = tile_names()
    assert tiles[-1] == "S36YY"


def test_tile_names_length():
    tiles = tile_names()
    assert len(tiles) == 32


def test_tile_pairs_first():
    pairs = tile_pairs(tile_names)
    assert pairs[0] == ("rf0XX", "S06XX")


def test_tile_pairs_last():
    pairs = tile_pairs(tile_names)
    assert pairs[-1] == ("rf1YY", "S36YY")


def test_tile_pairs_length():
    pairs = tile_pairs(tile_names)
    assert len(pairs) == 56


def test_time_tree_dates():
    dates, _ = time_tree("2020-01-01", "2020-01-02")
    assert dates == ["2020-01-01", "2020-01-02"]


def test_time_tree_time_stamps_start():
    _, time_stamps = time_tree("2020-01-01", "2020-01-02")
    assert time_stamps[0][0] == "2020-01-01-00:00"


def test_time_tree_time_stamps_end():
    _, time_stamps = time_tree("2020-01-01", "2020-01-02")
    assert time_stamps[-1][-1] == "2020-01-02-23:30"


def test_time_tree_time_stamps_size_day():
    _, time_stamps = time_tree("2020-01-01", "2020-01-02")
    assert len(time_stamps[0]) == 48
