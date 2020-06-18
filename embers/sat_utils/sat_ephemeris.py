import numpy as np
import skyfield as sf
from astropy.time import Time
from skyfield.api import Topos, load


# Skyfield Timescale
try:
    ts = load.timescale()
except:
    ts = load.timescale(builtin=True)

# Position of MWA site in Lat/Lon/Elevation
MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)


def load_tle(tle_path):
    """Loads a tle.
    
    Opens TLE file and uses skyfield to load each pair of TLE lines
    into an EarthSatellite object. Determines epoch for each.

    Agrs:
        the_path: Path to tle file

    Returns:
        sats: List of skyfield earth satellite objects
        epochs: Time at which each set of TLE lines is most accurate
    """

    # open a Two Line Element (TLE) file
    with open(tle_path, "r") as f:
        tle_list = [line.strip() for line in f.read().split("\n") if line is not ""]

    # List of satellite ephemeris from each set of TLE in the opened file
    # Epoch: Time at which TLE is most accurate (in JD - Julian Date)
    sats = []
    epochs = []

    for i in range(0, len(tle_list), 2):
        sat = sf.sgp4lib.EarthSatellite(tle_list[i], tle_list[i + 1])
        epoch = sat.model.jdsatepoch
        sats.append(sat)
        epochs.append(epoch)

    sats = np.asarray(sats)
    epochs = np.asarray(epochs)

    return (sats, epochs)


def epoch_ranges(epochs):
    """Creates a time array
    
    Time array with intervals corresponding to epochs of 
    best accuracy from the TLE pair of lines in TLE file.
    This is done to ensure that the most relevanrt TLE is
    used in the analysis.
     
     Args:
        epochs: List of epochs derived from TLE file

     Returns:
        epoch_range: List of times, between which TLE is most accurate
    """

    # Midpoints between epochs of successive TLEs
    midpoints = list((epochs[1:] + epochs[:-1]) / 2)
    epoch_range = [epochs[0]] + midpoints + [epochs[-1]]

    return epoch_range


def epoch_time_array(epoch_range, index_epoch, cadence):
    """Create a time array.
    
    Skyfield time object [time vector] at which
    the sat position will be determined.

    Args:
        index_epoch: Index of epoch range to be converted to time array
        t_step: Cadence of time array
        epoch_range: List of times, between which TLE is most accurate

    Returns:
        t_arr: Skyfield time object
        index_epoch: Same as above
    """

    # find time between epochs/ midpoints in seconds, using Astropy Time
    t1 = Time(epoch_range[index_epoch], scale="tt", format="jd").gps
    t2 = Time(epoch_range[index_epoch + 1], scale="tt", format="jd").gps
    dt = round(t2 - t1)

    t3 = Time(epoch_range[index_epoch], scale="tt", format="jd").iso
    date, time = t3.split()
    year, month, day = date.split("-")
    hour, minute, _ = time.split(":")

    # Create a time array, every 2[cadence] seconds, starting at the first epoch, approximately
    seconds = np.arange(0, dt, cadence)
    t_arr = ts.tt(int(year), int(month), int(day), int(hour), int(minute), seconds)

    return (t_arr, index_epoch)


def sat_pass(sats, t_arr, index_epoch):
    """Finds satellite passes.
    
    Calculates alt/az of sat at every position instant of time 
    array (t_arr). Finds all the times that the sat is above
    the horizon, and returns the index t_arr at which the sat
    rose and set.
    
    Args:
        sats: list of Skyfield EarthSatellite objects
        t_arr: Skyfield time object
        index_epoch: Position in epoch_range
    
    Returns:
        passes: 2D array with pairs of indicies corresponding to rise/set of sat
        alt: Array of altitudes of sat at t_arr times
        az: Array of azimuths of sat at t_arr times
    """

    # Define Satellite to be the first one in sats list
    # Find possition of sat at each timestep of time array

    if len(t_arr) > 0:

        satellite = sats[index_epoch]
        orbit = (satellite - MWA).at(t_arr)
        alt, az, _ = orbit.altaz()

        # Check if sat is above the horizon (above -1 degrees), return boolean array
        above_horizon = alt.degrees >= -1

        # Indicies of rare times that sats are above the horizon
        (indicies,) = above_horizon.nonzero()

        # Boundary times at which the sat either rises or sets
        (boundaries,) = np.diff(above_horizon).nonzero()

        if above_horizon[0] == True:
            boundaries = [indicies[0]] + list(boundaries)
            boundaries = np.asarray(boundaries)

        if above_horizon[-1] == True:
            boundaries = list(boundaries) + [indicies[-1]]
            boundaries = np.asarray(boundaries)

        if above_horizon[-1] == True and above_horizon[0] == True:
            boundaries = [indicies[0]] + list(boundaries) + [indicies[-1]]
            boundaries = np.asarray(boundaries)

        # Reshape into pairs rise & set indicies
        passes = boundaries.reshape(len(boundaries) // 2, 2)

        return (passes, alt, az)


def ephem_data(t_arr, pass_index, alt, az):
    """Satellite Ephemeris Data.
    
    creates rise time, set times and alt, az arrays.
    
    Args:
        t_arr: Skyfield time array object
        pass_index: One pair of sat indicies from passed 2D array
        alt: Array of altitudes of sat at t_arr times
        az: Array of azimuths of sat at t_arr times

    Returns:
        t_rise: Rise time of sat in gps seconds (Astropy)
        t_set: Set time of sat in gps seconds (Astropy)
        sat_alt: Altitude of satellite while it is above the horizon. Array.
        sat_az: Azimuth of satellite while it is above the horizon. Array.
    """

    i, j = pass_index

    # A list of times at which alt/az were calculated
    # Convert to unix time to match the rf explorer timestamps
    time_array = Time(t_arr[i : j + 1].tt, scale="tt", format="jd").unix

    theta = list(az.radians)
    r = list(alt.degrees)

    sat_az = theta[i : j + 1]
    sat_alt = r[i : j + 1]

    return (time_array, sat_alt, sat_az)


def sat_plot(sat_id, alt, az, num_passes, alpha=0.3):
    """Plots satellite passes
    
    Args:
        alt: list of altitude values
        az: list of azimuth values
        num_passes: Number of satellite passes
    """

    import matplotlib

    # Force matplotlib to not use X-Server backend
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Set up the polar plot.
    plt.style.use("dark_background")
    figure = plt.figure(figsize=(6, 6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0, 30, 60, 90], angle=22)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_title(f"Satellite {sat_id} Sky Coverage: {num_passes} Passes", y=1.08)
    ax.grid(color="#8bbabb", linewidth=1.6, alpha=0.6)
    plt.tight_layout()

    for i in range(len(alt)):
        plt.plot(az[i], alt[i], "-", linewidth=1.6, alpha=alpha, color="#1f4e5f")

    # Return plot for saving, showing, etc
    return plt


if __name__ == "__main__":

    import os
    import json
    import argparse

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
        default="./../../outputs/sat_ephemeris/TLE",
        help="Directory where TLE files are saved. Default=./../../outputs/sat_ephemeris/TLE",
    )
    parser.add_argument(
        "--cadence",
        metavar="\b",
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
    sat_name = args.sat
    tle_dir = args.tle_dir
    cadence = int(args.cadence)
    out_dir = args.out_dir

    tle_path = f"{tle_dir}/{sat_name}.txt"

    os.makedirs(os.path.dirname(out_dir), exist_ok=True)

    sat_ephem = {}
    sat_ephem["sat_id"] = [sat_name]
    sat_ephem["time_array"] = []
    sat_ephem["sat_alt"] = []
    sat_ephem["sat_az"] = []

    sats, epochs = load_tle(tle_path)
    epoch_range = epoch_ranges(epochs)

    for i in range(len(epoch_range) - 1):
        try:
            t_arr, index_epoch = epoch_time_array(epoch_range, i, cadence)
            passes, alt, az = sat_pass(sats, t_arr, index_epoch)

            for pass_index in passes:
                time_array, sat_alt, sat_az = ephem_data(t_arr, pass_index, alt, az)

                time_array = list(time_array)

                # We don't care about sat passes shorter than a minute (3*20sec)
                if len(time_array) >= 3:
                    sat_ephem["time_array"].append(time_array)
                    sat_ephem["sat_alt"].append(sat_alt)
                    sat_ephem["sat_az"].append(sat_az)
        except Exception:
            pass

    with open(f"{out_dir}/{sat_name}.json", "w") as outfile:
        json.dump(sat_ephem, outfile, indent=4)
