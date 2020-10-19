"""
Waterfall Batch
===============

Create waterfall plots for all rf data files recorded within a date interval.
Output files are saved to ``./embers_out/rf_tools/waterfalls``

"""

import argparse
import concurrent.futures
import logging
from itertools import repeat
from pathlib import Path

import pkg_resources
from embers.rf_tools.rf_data import batch_waterfall, tile_names, time_tree

_parser = argparse.ArgumentParser(
    description="""
    Create waterfall plots for all rf_files within a date interval
    """
)

_parser.add_argument(
    "--start_date", metavar="\b", default="", help="start date in YYYY-MM-DD format"
)

_parser.add_argument(
    "--stop_date", metavar="\b", default="", help="stop date in YYYY-MM-DD format"
)

_parser.add_argument(
    "--data_dir", metavar="\b", default="", help="root of dir where rf data is saved"
)

_parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/rf_tools",
    help="Dir where colormap sample plot is saved. Default=./embers_out/rf_tools",
)

_args = _parser.parse_args()
_start_date = _args.start_date
_stop_date = _args.stop_date
_data_dir = _args.data_dir
_out_dir = _args.out_dir

# if no input file provided, use sample package data
if _data_dir == "":
    print("----------------------------------------------------------")
    print("No data dir provided, using packaged sample data")
    _data_dir = pkg_resources.resource_filename("embers.kindle", "data/rf_data/")

if _start_date == "":
    print("No start_date provided, using 2019-10-01 for sample data")
    _start_date = "2019-10-01"

if _stop_date == "":
    print("No stop_date provided, using 2019-10-10 for sample data")
    print("")
    print(">>> waterfall_batch --help, for more options")
    print("----------------------------------------------------------")
    _stop_date = "2019-10-10"


# Logging config
_log_dir = Path(f"{_out_dir}/waterfalls")
_log_dir.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=f"{_out_dir}/waterfalls/waterfall_batch.log",
    level=logging.INFO,
    format="%(levelname)s: %(funcName)s: %(message)s",
)


def waterfall_batch(start_date, stop_date, data_dir, out_dir):
    """
    Save a series of waterfall plots in parallel.

    Parameters
    ----------
    :param start_date: date in style YYYY-MM-DD
    :type start_date: str
    :param stop_date: date in style YYYY-MM-DD
    :type stop_date: str
    :param data_dir: path to root of rf data dir
    :type data_dir: str
    :param out_dir: path to output dir
    :type out_dir: str

    """

    dates, time_stamps = time_tree(start_date, stop_date)

    for tile in tile_names():
        for day in range(len(dates)):

            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(
                    batch_waterfall,
                    repeat(tile),
                    time_stamps[day],
                    repeat(data_dir),
                    repeat(out_dir),
                )

            for result in results:
                logging.info(result)


def main():
    """Execute waterfall from terminal."""

    print(f"Processing rf data files between {_start_date} and {_stop_date}")
    print(f"Saving waterfall plots to: ./{_log_dir}")
    waterfall_batch(_start_date, _stop_date, _data_dir, _out_dir)
