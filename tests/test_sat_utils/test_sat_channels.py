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
