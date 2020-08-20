import shutil
from os import path
from pathlib import Path

from embers.rf_tools.align_data import (plot_savgol_interp, save_aligned,
                                        savgol_interp)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))

(ref_ali, tile_ali, time_array, _, _, _, _) = savgol_interp(
    f"{test_data}/rf_tools/rf_data/rf0XX/2019-10-01/rf0XX_2019-10-01-14:30.txt",
    f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt",
    savgol_window_1=11,
    savgol_window_2=15,
    polyorder=2,
    interp_type="cubic",
    interp_freq=1,
)

out_str = save_aligned(
    ("rf0XX", "S06XX"),
    "2019-10-01-14:30",
    11,
    15,
    2,
    "cubic",
    1,
    f"{test_data}/rf_tools/rf_data",
    f"{test_data}/rf_tools",
)


def test_savgol_interp_ref_tile_shape():
    assert ref_ali.shape == tile_ali.shape


def test_savgol_interp_ref_time_shape():
    assert time_array.shape[0] == 1779


def test_savgol_interp_ref_time_arrow():
    assert time_array[0] <= time_array[-1]


def test_plot_savgol_interp():
    plot_savgol_interp(
        ref=f"{test_data}/rf_tools/rf_data/rf0XX/2019-10-01/rf0XX_2019-10-01-14:30.txt",
        tile=f"{test_data}/rf_tools/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt",
        savgol_window_1=11,
        savgol_window_2=15,
        polyorder=2,
        interp_type="cubic",
        interp_freq=1,
        channel=14,
        out_dir=".",
    )
    png = Path("savgol_interp_sample.png")
    assert png.is_file() is True
    if png.is_file() is True:
        png.unlink()


def test_save_aligned_str():
    assert (
        out_str
        == f"Saved aligned file to {test_data}/rf_tools/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )


def test_save_aligned_file():
    ali_file = Path(
        f"{test_data}/rf_tools/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    assert ali_file.is_file() is True
    if ali_file.is_file() is True:
        shutil.rmtree(f"{test_data}/rf_tools/2019-10-01")


def test_save_aligned_err():
    out_str = save_aligned(
        ("rf0XX", "S06XX"),
        "2019-10-01-14:00",
        11,
        15,
        2,
        "cubic",
        1,
        f"{test_data}/rf_tools/rf_data",
        f"{test_data}/rf_tools",
    )
    assert type(out_str).__name__ == "FileNotFoundError"
