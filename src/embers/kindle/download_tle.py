"""
Download TLE
============

Download tle files for all satellies in :func:`~embers.sat_utils.sat_list.norad_ids`
within a date interval using Space-Track.org. TLE files are saved to
``./embers_out/sat_utils/TLE/*``

"""

import argparse
import os

from embers.sat_utils.sat_list import download_tle, norad_ids


def main():
    """
    Download tle files with the :func:`~embers.sat_utils.sat_list.download_tle` from space-tracks.org

    .. code-block:: console

        $ download_tle --help

    """
    _parser = argparse.ArgumentParser(
        description="""
        Downloads satellite TLE files from Space-tracks.org
        """
    )

    _parser.add_argument(
        "--start_date", metavar="\b", default="", help="start date in YYYY-MM-DD format"
    )

    _parser.add_argument(
        "--stop_date", metavar="\b", default="", help="stop date in YYYY-MM-DD format"
    )

    _parser.add_argument(
        "--st_ident", metavar="\b", default="", help="Space-Track.org login identity"
    )

    _parser.add_argument(
        "--st_pass", metavar="\b", default="", help="Space-Track.org login password"
    )

    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/sat_utils/TLE",
        help="Dir where satellite TLE files are saved. Default=./embers_out/sat_utils/TLE",
    )

    _args = _parser.parse_args()
    _start_date = _args.start_date
    _stop_date = _args.stop_date
    _st_ident = _args.st_ident
    _st_pass = _args.st_pass
    _out_dir = _args.out_dir

    n_ids = norad_ids()

    if _st_pass == "":
        # Check space-tracks.org credentials are saved as environment variables
        print(">>> download_tle --help, for usage details")
        try:
            _st_ident = os.environ.get("ST_USER")
            _st_pass = os.environ.get("ST_PASS")
        except Exception:
            pass

    print(f"TLE files saved to {_out_dir}")
    download_tle(
        _start_date,
        _stop_date,
        n_ids,
        st_ident=_st_ident,
        st_pass=_st_pass,
        out_dir=_out_dir,
    )
