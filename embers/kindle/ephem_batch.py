import json
import logging
import argparse
import numpy as np
import pkg_resources
from pathlib import Path
import concurrent.futures
from embers.sat_utils.sat_ephemeris import save_ephem

_parser = argparse.ArgumentParser(
    description="""
        Code which converts the TLE files downloaded with download_TLE.py
        into satellite ephemeris data: rise time, set time, alt/az arrays
        at a given time cadence. This is saved to a json file which will 
        be used to plot the satellite passes.
        """
)

_parser.add_argument(
    "--tle_dir", metavar="\b", default="", help="Path to directory with TLE files"
)
_parser.add_argument(
    "--cadence",
    metavar="\b",
    type=int,
    default=4,
    help="Rate at which sat alt/az is computed. default=4s",
)
_parser.add_argument(
    "--location",
    metavar="\b",
    type=json.loads,
    default=(-26.703319, 116.670815, 337.83),
    help="Geographic location where satellite ephemeris is to be determined. Default=MWA:(-26.703319, 116.670815, 337.83)",
)
_parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/sat_utils/",
    help="Path to output directory. Default=./embers_out/sat_utils/",
)

_args = _parser.parse_args()
_tle_dir = _args.tle_dir
_cadence = _args.cadence
_location = _args.location
_out_dir = _args.out_dir

# if no input file provided, use sample package data
if _tle_dir == "":
    print("----------------------------------------------------")
    print("No tle_dir path provided, using packaged sample data")
    print(">>> ephem_batch --help, for more options")
    print("----------------------------------------------------")
    _tle_file = Path(
        pkg_resources.resource_filename("embers.kindle", "data/TLE/25984.txt")
    )
    _tle_dir = _tle_file.parents[0]

# Logging config
_log_dir = Path(f"{_out_dir}/ephem_data")
_log_dir.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=f"{_out_dir}/ephem_data/ephem_batch.log",
    level=logging.INFO,
    format="%(levelname)s: %(funcName)s: %(message)s",
)


def ephem_batch(sat_name, tle_dir, cadence, location, out_dir):
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





def main():
    """Execute save ephem from terminal."""

    save_ephem(
        _sat_name,
        tle_dir=_tle_dir,
        cadence=_cadence,
        location=_location,
        out_dir=_out_dir,
    )
