"""
Align Batch
===========

Interpolate and Savitzky–Golay smooth a pair of raw rf data files
before temorally aligning the. Save aligned data as npz files to
``./embers_out/rf_tools/align_data``

"""

import argparse
import concurrent.futures
import logging
from itertools import repeat
from pathlib import Path

import pkg_resources
from embers.rf_tools.align_data import save_aligned
from embers.rf_tools.rf_data import tile_names, tile_pairs, time_tree

_parser = argparse.ArgumentParser(
    description="""
    Batch align pairs of rf data files between a date interval
    """
)

_parser.add_argument(
    "--start_date", metavar="\b", default="", help="start date in YYYY-MM-DD format"
)

_parser.add_argument(
    "--stop_date", metavar="\b", default="", help="stop date in YYYY-MM-DD format"
)

_parser.add_argument(
    "--savgol_window_1",
    metavar="\b",
    type=int,
    default=11,
    help="First savgol window. Default=11",
)

_parser.add_argument(
    "--savgol_window_2",
    metavar="\b",
    type=int,
    default=15,
    help="Second savgol window. Default=15",
)

_parser.add_argument(
    "--polyorder", metavar="\b", type=int, default=2, help="Polynomial order. Default=2"
)

_parser.add_argument(
    "--interp_type",
    metavar="\b",
    default="cubic",
    help="Interpolation type. Default=cubic",
)

_parser.add_argument(
    "--interp_freq",
    metavar="\b",
    type=int,
    default=1,
    help="Interpolation frequency. Default=1",
)

_parser.add_argument(
    "--data_dir", metavar="\b", default="", help="root of dir where rf data is saved"
)

_parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/rf_tools/align_data",
    help="Dir where savgol-interp sample plot is saved. Default=./embers_out/rf_tools/align_data",
)

_args = _parser.parse_args()
_start_date = _args.start_date
_stop_date = _args.stop_date
_savgol_window_1 = _args.savgol_window_1
_savgol_window_2 = _args.savgol_window_2
_polyorder = _args.polyorder
_interp_type = _args.interp_type
_interp_freq = _args.interp_freq
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
    print(">>> align_batch --help, for more options")
    print("----------------------------------------------------------")
    _stop_date = "2019-10-10"


# Logging config
_log_dir = Path(f"{_out_dir}")
_log_dir.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=f"{_out_dir}/align_batch.log",
    level=logging.INFO,
    format="%(levelname)s: %(funcName)s: %(message)s",
)


def align_batch(
    start_date=None,
    stop_date=None,
    savgol_window_1=None,
    savgol_window_2=None,
    polyorder=None,
    interp_type=None,
    interp_freq=None,
    data_dir=None,
    out_dir=None,
):

    dates, time_stamps = time_tree(start_date, stop_date)

    for pair in tile_pairs(tile_names()):

        for day in range(len(dates)):

            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(
                    save_aligned,
                    repeat(pair),
                    time_stamps[day],
                    repeat(savgol_window_1),
                    repeat(savgol_window_2),
                    repeat(polyorder),
                    repeat(interp_type),
                    repeat(interp_freq),
                    repeat(data_dir),
                    repeat(out_dir),
                )

            for result in results:
                logging.info(result)


def main():
    """Execute align_batch from terminal."""

    print(f"Aligned files saved to: {_out_dir}")
    align_batch(
        start_date=_start_date,
        stop_date=_stop_date,
        savgol_window_1=_savgol_window_1,
        savgol_window_2=_savgol_window_2,
        polyorder=_polyorder,
        interp_type=_interp_type,
        interp_freq=_interp_freq,
        data_dir=_data_dir,
        out_dir=_out_dir,
    )
