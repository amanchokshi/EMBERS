"""
MWA Pointings
-------------

Download MWA metadata and determine the pointings of rf observations

"""

import argparse

from embers.mwa_utils.mwa_pointings import mwa_point_meta


def main():
    """
    Download MWA pointing metadata using the :func:`~embers.mwa_utils.mwa_pointings.mwa_point_meta` function.

    .. code-block:: console

        $ mwa_pointings --help

    """

    _parser = argparse.ArgumentParser(
        description="""
        Download MWA metadata and determine the pointing of rf observations
        """
    )

    _parser.add_argument(
        "--start_date",
        metavar="\b",
        default="2019-10-09",
        help="start date in YYYY-MM-DD format. Default=2019-10-10",
    )
    _parser.add_argument(
        "--stop_date",
        metavar="\b",
        default="2019-10-11",
        help="stop date in YYYY-MM-DD format. Default=2019-10-10",
    )
    _parser.add_argument(
        "--num_pages",
        metavar="\b",
        default=4,
        type=int,
        help="Number of pages of metadata to download. Visit ws.mwatelescope.org/metadata/find. Default=2",
    )
    _parser.add_argument(
        "--time_thresh",
        metavar="\b",
        type=int,
        default=1,
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
        default="./tiles_data",
        help="Path to root of directory with raw rf data. Default=./tiles_data",
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

    mwa_point_meta(
        _start_date, _stop_date, _num_pages, _time_thresh, _time_zone, _rf_dir, _out_dir
    )
    print(f"MWA tile pointing data saved to {_out_dir}")
