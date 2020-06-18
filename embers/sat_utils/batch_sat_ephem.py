import os
import sys
import json
import argparse
from itertools import repeat
import sat_ephemeris as se
import sat_ids

import numpy as np
import skyfield as sf
from pathlib import Path
import concurrent.futures
from astropy.time import Time
from skyfield.api import Topos, load


# Save ephem date for one satellite
def sat_json(tle_dir, cadence, out_dir, sat_id):
    try:
        tle_path = f"{tle_dir}/{sat_id}.txt"

        sat_ephem = {}
        sat_ephem["sat_id"] = [sat_id]
        sat_ephem["time_array"] = []
        sat_ephem["sat_alt"] = []
        sat_ephem["sat_az"] = []

        sats, epochs = se.load_tle(tle_path)
        epoch_range = se.epoch_ranges(epochs)

        for i in range(len(epoch_range) - 1):
            try:
                t_arr, index_epoch = se.epoch_time_array(epoch_range, i, cadence)
                passes, alt, az = se.sat_pass(sats, t_arr, index_epoch)

                for pass_index in passes:
                    time_array, sat_alt, sat_az = se.ephem_data(
                        t_arr, pass_index, alt, az
                    )

                    time_array = list(time_array)

                    # We don't care about sat passes shorter than a minute (3*20sec)
                    if len(time_array) >= 3:
                        sat_ephem["time_array"].append(time_array)
                        sat_ephem["sat_alt"].append(sat_alt)
                        sat_ephem["sat_az"].append(sat_az)
            except Exception:
                pass
        with open(f"{out_dir}/{sat_id}.json", "w") as outfile:
            json.dump(sat_ephem, outfile, indent=4)

        return f"Saved {sat_id}.json"
    except Exception:
        return f"ERROR! Couldn't save {sat_id}.json. TLE was empty!"


if __name__ == "__main__":

    # Skyfield Timescale
    ts = load.timescale(builtin=True)

    # Position of MWA site in Lat/Lon/Elevation
    MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)

    parser = argparse.ArgumentParser(
        description="""
            Code which converts the TLE files downloaded with download_TLE.py
            into satellite ephemeris data: rise time, set time, alt/az arrays
            at a given time cadence. This is saved to a json file which will 
            be used to plot the satellite passes.
            """
    )

    parser.add_argument(
        "--tle_dir",
        metavar="\b",
        default="./../../outputs/sat_ephemeris/TLE",
        help="Directory where TLE files are saved. Default=./../../outputs/sat_ephemeris/TLE",
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
        default="./../../outputs/sat_ephemeris/ephem_json/",
        help="Path to output directory. Default=./../../outputs/sat_ephemeris/ephem_json/",
    )

    args = parser.parse_args()

    tle_dir = args.tle_dir
    cadence = int(args.cadence)
    out_dir = args.out_dir

    # Load list of satellites (Nodar ids)
    sat_list = list(sat_ids.norad_ids.values())

    # Save logs
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f"{out_dir}/sat_ephem_logs.txt", "a")

    # Parellization Magic Here!
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(
            sat_json, repeat(tle_dir), repeat(cadence), repeat(out_dir), sat_list
        )

    for result in results:
        print(result)
