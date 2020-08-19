import shutil
from os import path
from pathlib import Path

from embers.sat_utils.chrono_ephem import (interp_ephem, obs_times,
                                           save_chrono_ephem, write_json)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


def test_obs_times_human():
    obs_time, obs_unix, obs_unix_end = obs_times(
        "Australia/Perth", "2020-01-01", "2020-01-02"
    )
    assert obs_time[0] == "2020-01-01-00:00"


def test_obs_times_interval():
    obs_time, obs_unix, obs_unix_end = obs_times(
        "Australia/Perth", "2020-01-01", "2020-01-02"
    )
    assert obs_unix_end[0] - obs_unix[0] == 1800

#  def test_interp_ephem():
