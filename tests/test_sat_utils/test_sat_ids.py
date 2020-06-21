# Package imports
import pytest

# Embers imports
from embers.sat_utils.sat_ids import norad_ids


def test_norad_ids_length():
    assert len(norad_ids().values()) == 73


def test_norad_ids_first():
    assert norad_ids()["ORBCOMM-X"] == 21576


def test_norad_ids_last():
    assert norad_ids()["Meteor M2-2"] == 44387
