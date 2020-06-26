import json
import argparse
import numpy as np
from pathlib import Path
from embers.sat_utils.sat_ephemeris import save_ephem

parser = argparse.ArgumentParser(
    description="""
        Code which converts the TLE files downloaded with download_TLE.py
        into satellite ephemeris data: rise time, set time, alt/az arrays
        at a given time cadence. This is saved to a json file which will 
        be used to plot the satellite passes.
        """
)

parser.add_argument(
    "--sat", metavar="\b", help="The Norad cat ID of satellite. Example: 21576"
)
parser.add_argument(
    "--tle_dir",
    metavar="\b",
    default="../../embers_out/sat_utils/TLE/",
    help="Directory where TLE files are saved. Default=../../embers_out/sat_utils/TLE/",
)
parser.add_argument(
    "--cadence",
    metavar="\b",
    type=int,
    default=4,
    help="Rate at which sat alt/az is computed. Default=4s",
)
parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/sat_utils/",
    help="Path to output directory. Default=./embers_out/sat_utils/",
)

args = parser.parse_args()
sat_name = args.sat
tle_dir = args.tle_dir
cadence = args.cadence
out_dir = args.out_dir

# Position of MWA site in Lat/Lon/Elevation
MWA = (-26.703319, 116.670815, 337.83)

save_ephem(sat_name, tle_dir=tle_dir, cadence=cadence, location=MWA, out_dir=out_dir)

