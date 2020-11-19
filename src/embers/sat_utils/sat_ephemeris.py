"""
Satellite Ephemeris
-------------------

A set of tools to calculate satellite ephemeris
from TLE files.

"""

import concurrent.futures
import logging
from itertools import repeat
from pathlib import Path

import matplotlib as mpl
import numpy as np
import skyfield as sf
from astropy.time import Time
from matplotlib import pyplot as plt
from skyfield.api import Topos, load

mpl.use("Agg")


def load_tle(tle_file):
    """
    Extract orbital parameters from a TLE file.

    Instantiate an :class:`~skyfield.sgp4lib.EarthSatellite` for each pair of
    TLE lines in the TLE file, Also return the 'epoch' of each
    :class:`~skyfield.sgp4lib.EarthSatellite` object, which is the date and
    time for which the set of TLE lines is most accurate.

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
        tle_list = [line.strip() for line in f.read().split("\n") if line != ""]

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


def epoch_time_array(
    epoch_range, index_epoch=None, cadence=None,
):
    """Create a Skyfield :class:`~skyfield.timelib.Timescale` object at which to evaluate satellite positions.

    Begins by downloading up-to-date time files, using the skyfield :class:`~skyfield.iokit.Loader` class,
    needed to accurate converted between various time formats. See `Dates and Time
    <https://rhodesmill.org/skyfield/time.html>`_  for more info. The files are saved to
    :samp:`./embers_out/sat_utils/skyfield-data`.

    For a particular time inverval in :samp:`epoch_range`, chosen by :samp:`index_epoch`, a
    Skyfield :class:`~skyfield.timelib.Timescale` array object is generated, at a given time :samp:`cadence`. This time array will be used by :func:`~embers.sat_utils.sat_ephemeris.sat_pass` to compute the position of satellites at each time.

    .. code-block:: python

        from embers.sat_utils.sat_ephemeris import load_tle, epoch_ranges, epoch_time_array
        sats, epochs = load_tle('~/embers-data/TLE/21576.txt')
        epoch_range = epoch_ranges(epochs)
        index_epoch = 0     # select first time interval from epoch_range
        cadence = 10        # evaluate satellite position every 10 seconds

        t_arr, index_epoch = epoch_time_array(epoch_range, index_epoch, cadence)

    :param index_epoch: Index of :samp:`epoch_range` to be converted to time array :class:`~int`
    :param epoch_range: List of time intervals where an epoch is most accurate, from :func:`embers.sat_utils.sat_ephemeris.epoch_ranges`
    :param cadence: Time cadence at which to evaluate sat position, in seconds :class:`~int`

    :returns:
        A :class:`~tuple` of (t_arr, index_epoch)

        - t_arr: Skyfield :class:`~skyfield.timelib.Timescale` object with array of times at given :samp:`cadence`
        - index_epoch: Index of :samp:`epoch_range` to be converted to time array :class:`~int`

    """

    # Skyfield Timescale
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

    # Select satellite from sats with index_epoch
    # Find position of sat at each timestep of t_arr
    satellite = sats[index_epoch]
    orbit = (satellite - position).at(t_arr)
    alt, az, _ = orbit.altaz()

    if alt.degrees.shape[0] > 0:

        # Check if sat is above the horizon (above -1 degrees), return boolean array
        above_horizon = alt.degrees >= -1

        # Indicies of rare times that sats are above the horizon
        (indicies,) = above_horizon.nonzero()

        # Boundary times at which the sat either rises or sets
        (boundaries,) = np.diff(above_horizon).nonzero()

        if above_horizon[0]:
            boundaries = [indicies[0]] + list(boundaries)
            boundaries = np.asarray(boundaries)

        if above_horizon[-1]:
            boundaries = list(boundaries) + [indicies[-1]]
            boundaries = np.asarray(boundaries)

        # Reshape into pairs rise & set indicies
        passes = boundaries.reshape(len(boundaries) // 2, 2)

        return (passes, alt, az)

    else:
        return None


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
    time_array = Time(t_arr.tt[i : j + 1], scale="tt", format="jd").unix

    sat_az = az.radians[i : j + 1]
    sat_alt = alt.degrees[i : j + 1]

    return (time_array, sat_alt, sat_az)


def sat_plot(sat_id, alt, az, alpha=0.5):
    """Plots satellite passes

    :param sat_id: Norad catalogue ID :class:`~str`
    :param alt: :samp:`Altitude` :class:`~list`
    :param az: :samp:`Azimuth` :class:`~list`
    :param alpha: transparency of individual passes, default=0.5

    :return:
        - :func:`~matplotlib.pyplot.plot` object

    """

    plt.style.use("seaborn")
    figure = plt.figure(figsize=(6, 6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0, 30, 60, 90], angle=22)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.grid(color="#35635b", linewidth=1.8, alpha=0.7)

    for i in range(len(alt)):
        plt.plot(az[i], alt[i], "-", linewidth=1.8, alpha=alpha, color="#529471")

    ax.set_title(f"Satellite {sat_id} Sky Coverage: {len(alt)} Passes", y=1.08)
    plt.tight_layout()

    return plt


def save_ephem(sat, tle_dir, cadence, location, alpha, out_dir):
    """Save ephemeris of all satellite passes and plot sky coverage.

    This function brings everything in :mod:`~embers.sat_utils.sat_ephemeris` home.
    It converts a downloaded :samp:`TLE` file into arrays of :samp:`times`,
    :samp:`Altitudes` & :samp:`Azimuths` of when the satellite was above the
    horizon at a particular geographic :samp:`location`. These arrays are saved
    to the :samp:`out_dir` as an :class:`~numpy.savez_compressed` file. A plot
    of all satellite passes detected within the :samp:`TLE` file is also saved
    to the :samp:`out_dir`.

    .. code-block:: python

        from embers.sat_utils.sat_ephemeris import save_ephem
        sat="21576"
        cadence = 4
        tle_dir="~/embers-data/TLE"
        out_dir = "./embers_out/sat_utils"
        location = (-26.703319, 116.670815, 337.83) # MWA Telescope

        sat_ephem(sat, tle_dir, cadence=cadence, location, out_dir)

    .. code-block:: python

        Saved sky-coverage plot of sat [21576] to ./embers_out/sat_utils/ephem_plots
        Saved ephemeris of sat [21576] to ./embers_out/sat_utils/ephem_data

    :param sat: NORAD Catalogue ID of satellite :class:`~str`
    :param tle_dir: path to directory where :samp:`TLE` files are saved :class:`~str`
    :param cadence: time cadence at which to evaluate sat position, in seconds :class:`~int`
    :param location: The :samp:`gps` coordinates of the :samp:`location` at which satellite passes are to be computed. :samp:`location` is a :class:`~tuple` in the format (:samp:`latitude`, :samp:`longitude`, :samp:`elevation`), with :samp:`elevation` given in :samp:`meters`
    :param alpha: transparency of individual passes in :func:`~embers.sat_utils.sat_ephemeris.sat_plot` default=0.5
    :param out_dir: path to output directory :class:`~str`

    :returns:
        - satellite ephemeris at :samp:`location` and sky coverage ephemeris plot, saved to :samp:`out_dir`

    :raises FileNotFoundError: an input :samp:`TLE` file does not exist

    """

    # instantiate an empty dict
    sat_ephem = {}
    sat_ephem["sat_id"] = sat
    sat_ephem["time_array"] = []
    sat_ephem["sat_alt"] = []
    sat_ephem["sat_az"] = []

    # Make output directory tree
    Path(f"{out_dir}/ephem_data").mkdir(parents=True, exist_ok=True)
    Path(f"{out_dir}/ephem_plots").mkdir(parents=True, exist_ok=True)

    tle_path = Path(f"{tle_dir}/{sat}.txt")

    # Skip if the file is empty
    if tle_path.stat().st_size != 0:

        # Check if tle file exists
        tle_path.is_file()
        sats, epochs = load_tle(tle_path)
        epoch_range = epoch_ranges(epochs)

        for i in range(len(epoch_range) - 1):
            t_arr, index_epoch = epoch_time_array(
                epoch_range, index_epoch=i, cadence=cadence
            )

            try:
                passes, alt, az = sat_pass(sats, t_arr, index_epoch, location=location)

                for pass_index in passes:
                    time_array, sat_alt, sat_az = ephem_data(t_arr, pass_index, alt, az)

                    sat_ephem["time_array"].append(time_array)
                    sat_ephem["sat_alt"].append(sat_alt)
                    sat_ephem["sat_az"].append(sat_az)

            # Catch exceptions in sat_pass
            # sometimes sat_object is empty and can't be iterated over
            except Exception:
                pass

        plt = sat_plot(sat, sat_ephem["sat_alt"], sat_ephem["sat_az"], alpha=alpha)
        plt.savefig(f"{out_dir}/ephem_plots/{sat}.png")
        plt.close()
        np.savez_compressed(f"{out_dir}/ephem_data/{sat}.npz", **sat_ephem)

        return f"Saved sky coverage plot of satellite [{sat}] to {out_dir}ephem_plots/{sat}.png \nSaved ephemeris of satellite [{sat}] to {out_dir}ephem_data/{sat}.npz"

    return f"File {tle_dir}/{sat} is empty, skipping"


def ephem_batch(tle_dir, cadence, location, alpha, out_dir, max_cores=None):
    """
    Process ephemeris for multiple satellites in parallel.

    :param tle_dir: Path to directory where :samp:`TLE` files are saved :class:`~str`
    :param cadence: Time cadence at which to evaluate sat position, in seconds :class:`~int`
    :param location: The :samp:`gps` coordinates of the :samp:`location` at which satellite passes are to be computed. :samp:`location` is a :class:`~tuple` in the format (:samp:`latitude`, :samp:`longitude`, :samp:`elevation`), with :samp:`elevation` given in :samp:`meters`
    :param alpha: Transparency of individual passes in :func:`~embers.sat_utils.sat_ephemeris.sat_plot` default=0.5
    :param out_dir: Path to output directory :class:`~str`
    :param max_cores: Maximum number of cores to be used by this script. Default=None, which means that all available cores are used

    :returns:
        - satellite ephemeris at :samp:`location` and sky coverage ephemeris plot, saved to :samp:`out_dir`

    """

    # Logging config
    log_dir = Path(f"{out_dir}/ephem_data")
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=f"{out_dir}/ephem_data/ephem_batch.log",
        level=logging.INFO,
        format="%(levelname)s: %(funcName)s: %(message)s",
    )
    sat_names = [tle.stem for tle in Path(tle_dir).glob("*.txt")]

    with concurrent.futures.ProcessPoolExecutor(max_workers=max_cores) as executor:
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
