import json
import shutil
from os import path
from pathlib import Path

import numpy as np
from embers.sat_utils.sat_channels import (batch_window_map, good_chans,
                                           noise_floor, plt_channel, plt_sats,
                                           plt_window_chans, read_aligned,
                                           time_filter, window_chan_map)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))

ali_file = Path(
    f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
)

ali_file_2 = Path(
    f"{test_data}/rf_tools/align_data/2019-10-10/2019-10-10-02:30/rf0XX_S06XX_2019-10-10-02:30_aligned.npz"
)

chrono_file = Path(f"{test_data}/sat_utils/chrono_json/2019-10-01-14:30.json")
chrono_file_2 = Path(f"{test_data}/sat_utils/chrono_json/2019-10-10-02:30.json")


def test_read_aligned_ref_tile_shape():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    assert ref_pow.shape == tile_pow.shape


def test_read_aligned_times():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    assert times.shape[0] == 1779


def test_noise_floor():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    noise_threshold = noise_floor(1, 3, ref_pow)
    assert round(noise_threshold) == -104.0


def test_time_filter_I():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    intvl = time_filter(1569911400, 1569913100, times)
    assert intvl == [0, 1696]


def test_time_filter_II():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    intvl = time_filter(1569911410, 1569913100, times)
    assert intvl == [6, 1696]


def test_time_filter_III():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    intvl = time_filter(1569911405, 1569913200, times)
    assert intvl == [1, 1778]


def test_time_filter_IV():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    intvl = time_filter(1569911400, 1569913200, times)
    assert intvl == [0, 1778]


def test_time_filter_V():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    intvl = time_filter(1569911300, 1569911400, times)
    assert intvl is None


def test_plt_window_chans():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    plt = plt_window_chans(
        ref_pow, 12345, 1569911700, 1569912000, "Spectral", chs=[0, 4, 7, 14], good_ch=7
    )
    assert type(plt).__name__ == "module"


def test_plt_channel():
    ali_file = Path(
        f"{test_data}/rf_tools/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz"
    )
    ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
    noise_threshold = noise_floor(1, 3, ref_pow)
    plt = plt_channel(
        times, ref_pow[:, 46], np.median(ref_pow), 46, [-120, 5], noise_threshold, 30
    )
    assert type(plt).__name__ == "module"


def test_plt_sats():
    plt = plt_sats(
        ["25115", "25116"],
        f"{test_data}/sat_utils/chrono_json/2019-10-01-14:30.json",
        "2019-10-01-14:30",
    )
    assert type(plt).__name__ == "module"


def test_good_chans_41():
    good_chan = good_chans(
        ali_file,
        chrono_file,
        "41184",
        1,
        3,
        15,
        0.8,
        "2019-10-01-14:30",
        f"{test_data}/sat_utils/good_chans_tmp",
    )
    assert good_chan == 41


def test_good_chans_41_plts():
    good_chans(
        ali_file,
        chrono_file,
        "41184",
        1,
        3,
        15,
        0.8,
        "2019-10-01-14:30",
        f"{test_data}/sat_utils/good_chans_tmp",
        plots=True,
    )
    waterfall = Path(
        f"{test_data}/sat_utils/good_chans_tmp/window_plots/2019-10-01/2019-10-01-14:30/41184_waterfall_41.png"
    )

    assert waterfall.is_file()
    shutil.rmtree(f"{test_data}/sat_utils/good_chans_tmp")


def test_good_chans_45():
    good_chan = good_chans(
        ali_file_2,
        chrono_file_2,
        "41188",
        1,
        3,
        15,
        0.8,
        "2019-10-10-02:30",
        f"{test_data}/sat_utils/good_chans_tmp",
    )
    assert good_chan == 45


def test_good_chans_45_plts():
    good_chans(
        ali_file_2,
        chrono_file_2,
        "41188",
        1,
        3,
        15,
        0.8,
        "2019-10-10-02:30",
        f"{test_data}/sat_utils/good_chans_tmp",
        plots=True,
    )
    waterfall = Path(
        f"{test_data}/sat_utils/good_chans_tmp/window_plots/2019-10-10/2019-10-10-02:30/41188_waterfall_45.png"
    )

    assert waterfall.is_file()
    shutil.rmtree(f"{test_data}/sat_utils/good_chans_tmp")


def test_good_chans_46():
    good_chan = good_chans(
        ali_file,
        chrono_file,
        "25417",
        1,
        3,
        15,
        0.8,
        "2019-10-01-14:30",
        f"{test_data}/sat_utils/good_chans_tmp",
    )
    assert good_chan == 46


def test_good_chans_46_plts():
    good_chans(
        ali_file,
        chrono_file,
        "25417",
        1,
        3,
        15,
        0.8,
        "2019-10-01-14:30",
        f"{test_data}/sat_utils/good_chans_tmp",
        plots=True,
    )
    waterfall = Path(
        f"{test_data}/sat_utils/good_chans_tmp/window_plots/2019-10-01/2019-10-01-14:30/25417_waterfall_46.png"
    )

    assert waterfall.is_file()
    shutil.rmtree(f"{test_data}/sat_utils/good_chans_tmp")


def test_good_chans_plts_nochans():
    good_chans(
        ali_file,
        chrono_file,
        "44028",
        1,
        3,
        15,
        0.8,
        "2019-10-01-14:30",
        f"{test_data}/sat_utils/good_chans_tmp",
        plots=True,
    )
    waterfall = Path(
        f"{test_data}/sat_utils/good_chans_tmp/window_plots/2019-10-01/2019-10-01-14:30/44028_waterfall_window.png"
    )

    assert waterfall.is_file()
    shutil.rmtree(f"{test_data}/sat_utils/good_chans_tmp")


def test_window_chan_map():
    window_chan_map(
        f"{test_data}/rf_tools/align_data",
        f"{test_data}/sat_utils/chrono_json",
        1,
        3,
        15,
        0.8,
        "2019-10-01-14:30",
        f"{test_data}/sat_utils/good_chans_tmp",
        True,
    )
    chan_map = Path(
        f"{test_data}/sat_utils/good_chans_tmp/window_maps/2019-10-01-14:30.json"
    )

    assert chan_map.is_file()
    shutil.rmtree(f"{test_data}/sat_utils/good_chans_tmp")


def test_batch_window_map():
    batch_window_map(
        "2019-10-01",
        "2019-10-01",
        f"{test_data}/rf_tools/align_data",
        f"{test_data}/sat_utils/chrono_json",
        1,
        3,
        15,
        0.8,
        f"{test_data}/sat_utils/good_chans_tmp",
        plots=False,
    )
    chan_map = Path(
        f"{test_data}/sat_utils/good_chans_tmp/window_maps/2019-10-01-14:30.json"
    )
    assert chan_map.is_file()
    shutil.rmtree(f"{test_data}/sat_utils/good_chans_tmp")
