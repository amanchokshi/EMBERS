"""
MWA Pointings
-------------

Download MWA metadata and determine the pointings of rf observations

"""

import argparse

from embers.mwa_utils.mwa_pointings import mwa_point_meta

_parser = argparse.ArgumentParser(
    description="""
    Download MWA metadata and determine the pointing of rf observations
    """
)

_parser.add_argument(
    "--start_date", metavar="\b", default="", help="start date in YYYY-MM-DD format"
)
_parser.add_argument(
    "--stop_date", metavar="\b", default="", help="stop date in YYYY-MM-DD format"
)
_parser.add_argument(
    "--num_pages",
    metavar="\b",
    type=int,
    help="Number of pages of metadata to download. Visit ws.mwatelescope.org/metadata/find to find this",
)
_parser.add_argument(
    "--time_thresh",
    metavar="\b",
    type=int,
    default=100,
    help="Minimum integration in hours, at a pointing",
)
_parser.add_argument(
    "--time_zone",
    metavar="\b",
    default="Australia/Perth",
    help="pytz time zone where experiment was conducted. Default=Australia/Perth",
)
_parser.add_argument(
    "--rf_dir",
    metavar="\b",
    default="",
    help="Path to root of directory with raw rf data",
)
_parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/mwa_utils",
    help="Dir where MWA metadata will be saved. Default=./embers_out/mwa_utils",
)

_args = _parser.parse_args()
_start_date = _args.start_date
_stop_date = _args.stop_date
_num_pages = _args.num_pages
_time_thresh = _args.time_thresh
_time_zone = _args.time_zone
_rf_dir = _args.rf_dir
_out_dir = _args.out_dir


# if no input file provided, use sample package data
if _rf_dir == "":
    print("----------------------------------------------------------")
    print("required arguments missing")
    print(">>> mwa_pointings --help, for more options")
    print("----------------------------------------------------------")


def main():
    """Execute mwa_pointings from terminal."""

    mwa_point_meta(
        _start_date, _stop_date, _num_pages, _time_thresh, _time_zone, _rf_dir, _out_dir
    )
