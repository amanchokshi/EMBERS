import json
import logging
import argparse
import numpy as np
import pkg_resources
from pathlib import Path
import concurrent.futures
from itertools import repeat
from embers.sat_utils.sat_list import norad_ids
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
    "--alpha",
    metavar="\b",
    type=int,
    default=0.5,
    help="Alpha value for sky coverage plot. Defaut: 0.5. If too many satellite, reduce value",
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
_alpha = _args.alpha
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


def ephem_batch(tle_dir, cadence, location, alpha, out_dir):
    """
    Process ephemeris for multiple satellites in parallel.

    """
    sat_names = list(norad_ids().values())
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(
            save_ephem,
            sat_names,
            repeat(tle_dir),
            repeat(cadence),
            repeat(location),
            repeat(alpha),
            repeat(out_dir),
        )

    for result in results:
        logging.info(result)


def main():
    """Execute batch save_ephem from terminal."""

    print(f"Saving logs to {_out_dir}ephem_data")
    print(f"Saving sky coverage plots to {_out_dir}ephem_plots")
    print(f"Saving ephemeris of satellites to {_out_dir}ephem_data")
    ephem_batch(_tle_dir, _cadence, _location, _alpha, _out_dir)
