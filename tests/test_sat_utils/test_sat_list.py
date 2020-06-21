# Package imports
import pytest

# Embers imports
from embers.sat_utils.sat_list import norad_ids, download_tle


def test_norad_ids_length():
    assert len(norad_ids().values()) == 73


def test_norad_ids_first():
    assert norad_ids()["ORBCOMM-X"] == 21576


def test_norad_ids_last():
    assert norad_ids()["Meteor M2-2"] == 44387


def test_download_tle_fail(capfd):
    n_ids = norad_ids()
    download_tle("2020-01-01", "2020-02-01", n_ids, st_ident=None, st_pass=None, out_dir="./")
    captured = capfd.readouterr()
    assert captured.out == "Space-Track.org credentials not provided. Make an account before downloading TLEs\n"
