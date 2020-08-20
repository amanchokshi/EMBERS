from pathlib import Path

from embers.sat_utils.sat_list import download_tle, norad_ids


def test_norad_ids_length():
    assert len(norad_ids().values()) == 73


def test_norad_ids_first():
    assert norad_ids()["ORBCOMM-X"] == 21576


def test_norad_ids_last():
    assert norad_ids()["Meteor M2-2"] == 44387


def test_download_tle_fail(capfd):
    n_ids = norad_ids()
    download_tle(
        "2020-01-01", "2020-02-01", n_ids, st_ident=None, st_pass=None, out_dir="./"
    )
    captured = capfd.readouterr()
    assert (
        captured.out
        == "Space-Track.org credentials not provided. Make an account before downloading TLEs\n"
    )


def test_download_tle_http_err():
    n_ids = norad_ids()
    try:
        download_tle(
            "2020-01-01",
            "2020-02-01",
            n_ids,
            st_ident="test",
            st_pass="user",
            out_dir="./",
        )
    except Exception as e:
        assert type(e).__name__ == "HTTPError"


def test_download_tle_write():
    n_ids = {"ORBCOMM-X": 21576}
    print("Yellow")
    download_tle(
        "2020-01-01",
        "2020-02-01",
        n_ids,
        st_ident="test",
        st_pass="user",
        out_dir="./",
        sleep=0,
        mock=True,
    )
    txt = Path("21576.txt")
    assert txt.is_file() is True
    if txt.is_file() is True:
        txt.unlink()
