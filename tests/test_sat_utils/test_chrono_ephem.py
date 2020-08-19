import json
import shutil
from os import path
from pathlib import Path

import numpy as np
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


def test_interp_ephem_shape():
    ephem_data = Path(f"{test_data}/sat_utils/ephem_data/25986.npz")

    sat_ephem = np.load(ephem_data, allow_pickle=True)
    t_array = sat_ephem["time_array"]
    s_alt = sat_ephem["sat_alt"]
    s_az = sat_ephem["sat_az"]
    #  s_id = str(sat_ephem["sat_id"])

    time_interp, sat_alt, sat_az = interp_ephem(
        t_array[0], s_alt[0], s_az[0], "cubic", 1
    )
    assert len(time_interp) == 931


def test_interp_ephem_sky():
    ephem_data = Path(f"{test_data}/sat_utils/ephem_data/25986.npz")

    sat_ephem = np.load(ephem_data, allow_pickle=True)
    t_array = sat_ephem["time_array"]
    s_alt = sat_ephem["sat_alt"]
    s_az = sat_ephem["sat_az"]
    #  s_id = str(sat_ephem["sat_id"])

    time_interp, sat_alt, sat_az = interp_ephem(
        t_array[0], s_alt[0], s_az[0], "cubic", 1
    )
    assert len(sat_alt) == len(sat_az)


def test_write_json():
    data = [1, 2, 3]
    write_json(data, "tmp.json", f"{test_data}/sat_utils")
    tmp_path = Path(f"{test_data}/sat_utils/tmp.json")
    if tmp_path.is_file():
        with open(tmp_path, "r") as f:
            d_read = json.load(f)

            assert d_read == data
            tmp_path.unlink()


def test_save_chrono_ephem():
    out_dir = Path(f"{test_data}/sat_utils/ephem_chrono_tmp")
    save_chrono_ephem(
        "Australia/Perth",
        "2019-10-01",
        "2019-10-01",
        "cubic",
        1,
        f"{test_data}/sat_utils/ephem_data",
        out_dir,
    )
    files = [f.stem for f in out_dir.glob("*.json")]

    assert len(files) == 48
    shutil.rmtree(out_dir)
