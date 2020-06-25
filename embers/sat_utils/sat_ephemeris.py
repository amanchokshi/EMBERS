"""
Satellite Ephemeris
-------------------

A set of tools to calculate satellite ephemeris 
from TLE files. 

"""

import numpy as np
import skyfield as sf
from astropy.time import Time
import matplotlib.pyplot as plt
from skyfield.api import Topos, Loader


def load_tle(tle_file):
    """
    Extract orbital parameters from a TLE file.
    
    Instantiate an :class:`~skyfield.sgp4lib.EarthSatellite` for each pair of TLE lines in the TLE file, 
    Also return the 'epoch' of each :class:`~skyfield.sgp4lib.EarthSatellite` object, which is the date and time for which the set of TLE lines is most accurate.

    .. code-block:: python
        
        from embers.sat_utils.sat_ephemeris import load_tle
        sats, epochs = load_tle('~/embers-data/TLE/21576.txt')

    :param tle_file: path to TLE file :class:`~str`

    :returns:
        A :class:`~tuple` (sats, epochs)
        
        - sats - list of :class:`~skyfield.sgp4lib.EarthSatellite` objects, one for each pair of TLE lines
        - epochs - Julian date at which each set of TLE lines is most accurate

    """

    # open a Two Line Element (TLE) file
    with open(tle_file, "r") as f:
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
    """Optimise time intervals to make the most of different epochs of TLE pairs
    
    Creates a list with of times [epoch_range], pairs of successive elements
    of which correspond to time intervals at which a particular epoch has best 
    accuracy, before the next epoch becomes more accurate. This is done to 
    ensure that the most relevanrt TLE is used in the analysis.
    
    .. code-block:: python
        
        from embers.sat_utils.sat_ephemeris import load_tle, epoch_ranges
        sats, epochs = load_tle('~/embers-data/TLE/21576.txt')
        epoch_range = epoch_ranges(epochs)
     
    :param epochs: List of epochs from :func:`~embers.sat_utils.sat_ephemeris.load_tle`

    :return:
        epoch_range: List of times, between which each pair of TLEs is most accurate

    """

    # Midpoints between epochs of successive TLEs
    midpoints = list((epochs[1:] + epochs[:-1]) / 2)
    epoch_range = [epochs[0]] + midpoints + [epochs[-1]]

    return epoch_range


def epoch_time_array(epoch_range, index_epoch, cadence):
    """Create a Skyfield :class:`~skyfield.timelib.Timescale` object at which to evaluate satellite positions.
    
    Begins by downloading up-to-date time files, using the skyfield :class:`~skyfield.iokit.Loader` class, 
    needed to accurate converted between various time formats. See `Dates and Time 
    <https://rhodesmill.org/skyfield/time.html>`_  for more info. The files are saved to 
    :samp:`./embers_out/sat_utils/skyfield-data`. 

    For a particular time inverval in :samp:`epoch_range`, chosen by :samp:`index_epoch`, a 
    Skyfield :class:`~skyfield.timelib.Timescale` array object is generated, at a given time 
    :samp:`cadence`. This time array will be used by :func:`~embers.sat_utils.sat_ephemeris.sat_pass` 
    to compute the position of satellites at each time.  
    
    .. code-block:: python
        
        from embers.sat_utils.sat_ephemeris import load_tle, epoch_ranges, epoch_time_array
        sats, epochs = load_tle('~/embers-data/TLE/21576.txt')
        epoch_range = epoch_ranges(epochs)
        index_epoch = 0     # select first time interval from epoch_range
        cadence = 10        # evaluate satellite position every 10 seconds
        
        t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch, cadence)

    .. code-block:: console

        >>> skyfield time files downloaded to ./embers_out/sat_utils/skyfield-data
        >>> [#################################] 100% deltat.data
        >>> [#################################] 100% deltat.preds
        >>> [#################################] 100% Leap_Second.dat

    :param index_epoch: Index of :samp:`epoch_range` to be converted to time array :class:`~int`
    :param epoch_range: List of times intervals where an epoch is most accurate, from :func:`embers.sat_utils.sat_ephemeris.epoch_ranges` 
    :param cadence: time cadence at which to evaluate sat position, in seconds :class:`~int`

    :returns:
        A :class:`~tuple` of (t_arr, index_epoch)

        - t_arr: Skyfield :class:`~skyfield.timelib.Timescale` object with array of times at given :samp:`cadence`
        - index_epoch: Index of :samp:`epoch_range` to be converted to time array :class:`~int`

    """

    # Skyfield Timescale
    load = Loader("./embers_out/sat_utils/skyfield-data")

    try:
        ts = load.timescale()
    except Exception as e:
        if e != None:
            ts = load.timescale(builtin=True)

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


def sat_pass(sats, t_arr, index_epoch, location=None):
    """Find when a satellite passes above the horizon at a gps location.
    
    Calculate the :samp:`Altitude` & :samp:`Azimuth` of a 
    :class:`~skyfield.sgp4lib.EarthSatellite` object from 
    :func:`~embers.sat_utils.sat_ephemeris.load_tle` at 
    every instant of time in :samp:`t_arr` from  
    :func:`~embers.sat_utils.sat_ephemeris.epoch_time_array`. 
    Determine all the times that the satellite is above the 
    horizon, at a given gps :samp:`location` and returns the 
    pair of indices of :samp:`t_arr` at which the satellite rose 
    and set. 
    
    .. code-block:: python
        
        from embers.sat_utils.sat_ephemeris import load_tle, epoch_ranges, epoch_time_array, sat_pass
        sats, epochs = load_tle('~/embers-data/TLE/21576.txt')
        epoch_range = epoch_ranges(epochs)
        index_epoch = 0     # select first time interval from epoch_range
        cadence = 10        # evaluate satellite position every 10 seconds
        t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch, cadence)
        MWA = (-26.703319, 116.670815, 337.83)   # gps coordinates of MWA Telescope

        passes, alt, az = sat_pass(sats, t_arr, index_epoch, location=MWA)
    
    :param sats: list of :class:`~skyfield.sgp4lib.EarthSatellite` objects
    :param t_arr: skyfield :class:`~skyfield.timelib.Timescale` object with array of times
    :param index_epoch: Index of :samp:`epoch_range` :class:`~int`
    :param location: The :samp:`gps` coordinates of the :samp:`location` at which satellite passes are to be computed. :samp:`location` is a :class:`~tuple` in the format (:samp:`latitude`, :samp:`longitude`, :samp:`elevation`), with :samp:`elevation` given in :samp:`meters` 
    
    :returns:
        A :class:`~tuple` of (passes, alt, az)
        
        - passes: 2D array with pairs of indicies of :samp:`t_arr` corresponding to rise/set of satellite :class:`~numpy.ndarray`
        - alt: Array of :samp:`Altitudes` of sat at :samp:`t_arr` times :class:`~numpy.ndarray`
        - az: Array of :samp:`Azimuths` of sat at :samp:`t_arr` times :class:`~numpy.ndarray`

    """

    # Position where sat passes are to be determined in Lat/Lon/Elevation
    position = Topos(
        latitude=location[0], longitude=location[1], elevation_m=location[2]
    )

    if len(t_arr) > 0:

        # Select satellite from sats with index_epoch
        # Find position of sat at each timestep of t_arr
        satellite = sats[index_epoch]
        orbit = (satellite - position).at(t_arr)
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
    """Satellite Ephemeris data (time, alt, az arrays ) for a single satellite pass.
    
    .. code-block:: python
        
        from embers.sat_utils.sat_ephemeris import load_tle, epoch_ranges, epoch_time_array, sat_pass, ephem_data
        sats, epochs = load_tle('~/embers-data/TLE/21576.txt')
        epoch_range = epoch_ranges(epochs)
        index_epoch = 0     # select first time interval from epoch_range
        cadence = 10        # evaluate satellite position every 10 seconds
        t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch, cadence)
        MWA = (-26.703319, 116.670815, 337.83)   # gps coordinates of MWA Telescope
        passes, alt, az = sat_pass(sats, t_arr, index_epoch, location=MWA)

        time_array, sat_alt, sat_az = ephem_data(t_arr, passes[0], alt, az)
    
    :param t_arr: skyfield :class:`~skyfield.timelib.Timescale` object with array of times
    :param pass_index: One pair of sat indicies from :samp:`passes`
    :param alt: Array of altitudes of sat at :samp:`t_arr` times
    :param az: Array of azimuths of sat at :samp:`t_arr` times

    :returns:
        A :class:`~tuple` (time_array, sat_alt, sat_az)
        
        - time_array: times at which sat position is calculated :class:`~numpy.ndarray`
        - sat_alt: :samp:`Altitude` of satellite while it is above the horizon :class:`~numpy.ndarray`
        - sat_az: :samp:`Azimuth` of satellite while it is above the horizon :class:`~numpy.ndarray`
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


def sat_plot(sat_id, alt, az, num_passes, alpha=0.5):
    """Plots satellite passes
    
    :param sat_id: Norad catalogue ID :class:`~str`
    :param alt: :samp:`Altitude` :class:`~list`
    :param az: :samp:`Azimuth` :class:`~list`
    :param alpha: transparency of individual passes, default=0.5

    :return:
        - plt - :func:`~matplotlib.pyplot.plot` object
    
    """

    plt.style.use("seaborn")
    figure = plt.figure(figsize=(6, 6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0, 30, 60, 90], angle=22)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.grid(color="#8bbabb", linewidth=1.6, alpha=0.6)

    for i in range(len(alt)):
        plt.plot(az[i], alt[i], "-", linewidth=1.6, alpha=alpha, color="#1f4e5f")

    ax.set_title(f"Satellite {sat_id} Sky Coverage: {len(alt)} Passes", y=1.08)
    plt.tight_layout()

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

    # Position of MWA site in Lat/Lon/Elevation
    MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)

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
