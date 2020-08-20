import json
import shutil
from os import path
from pathlib import Path

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
ref_pow, tile_pow, times = read_aligned(ali_file=ali_file)
noise_threshold = noise_floor(1, 3, ref_pow)
print(round(noise_threshold))


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
