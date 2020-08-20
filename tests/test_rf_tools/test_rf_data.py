import shutil
from os import path
from pathlib import Path

from embers.rf_tools.rf_data import (batch_waterfall, plt_waterfall, read_data,
                                     single_waterfall, tile_names, tile_pairs,
                                     time_tree)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


def test_read_data_power_shape():
    power, _ = read_data(
        rf_file=f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert power.shape == (16655, 112)


def test_read_data_times_shape():
    _, times = read_data(
        rf_file=f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert list(times.shape)[0] == 16655


def test_read_data_power_times_shape():
    power, times = read_data(
        rf_file=f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert list(times.shape)[0] == list(power.shape)[0]


def test_read_data_return_type():
    _, times = read_data(
        rf_file=f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    assert type(times).__name__ == "ndarray"


def test_read_data_time_arrow():
    _, times = read_data(
        rf_file=f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
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


def test_plt_waterfall():
    power, times = read_data(
        rf_file=f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    plt = plt_waterfall(power, times, "S06XX-test")
    assert type(plt).__name__ == "module"


def test_plt_waterfall_savefig():
    power, times = read_data(
        rf_file=f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )
    plt = plt_waterfall(power, times, "S06XX-test")
    plt.savefig(f"{test_data}/rf_tools/plt_waterfall_test.png")
    test_file_path = Path(f"{test_data}/rf_tools/plt_waterfall_test.png")
    assert test_file_path.is_file() is True
    if test_file_path.is_file():
        test_file_path.unlink()


def test_single_waterfall():
    single_waterfall(
        f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt",
        f"{test_data}/rf_tools",
    )
    single_waterfall_png = Path(f"{test_data}/rf_tools/S06XX_2019-10-01-14:30.png")
    assert single_waterfall_png.is_file() is True
    if single_waterfall_png.is_file() is True:
        single_waterfall_png.unlink()


def test_batch_waterfall():
    batch_waterfall(
        "S06XX",
        "2019-10-01-14:30",
        f"{test_data}/rf_tools/rf_data",
        f"{test_data}/rf_tools",
    )
    batch_waterfall_png = Path(
        f"{test_data}/rf_tools/waterfalls/2019-10-01/2019-10-01-14:30/S06XX_2019-10-01-14:30.png"
    )
    assert batch_waterfall_png.is_file() is True
    if batch_waterfall_png.is_file() is True:
        shutil.rmtree(f"{test_data}/rf_tools/waterfalls")


def test_batch_waterfall_err():
    e = batch_waterfall("S06XX", "2019-10-01-14:30", ".", ".")
    assert type(e).__name__ == "FileNotFoundError"
