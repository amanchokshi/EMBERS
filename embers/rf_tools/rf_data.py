"""
RF Data
=======

A set of tools to decode raw rf data recored by
RF Explorers and visualise waterfall plots.

"""

import re
import time
import numpy as np
from pathlib import Path
from itertools import product
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from embers.rf_tools.colormaps import spectral

_spec, _ = spectral()


def read_data(rf_file=None):
    """
    Convert raw rf binary data into a power array and extracts a list of UNIX timestamps.

    Parameters
    ----------
    :param rf_file: path to a rf data file
    :type rf_file: str

    Returns
    -------
    :returns:
        - power - array of rf powers in dBm
        - times - array of UNIX times

    :rtype: (float, numpy.array(float))

    """

    with open(rf_file, "rb") as f:
        next(f)
        lines = f.readlines()

        times = []
        data_lines = []

        for line in lines:
            time, data = line.split("$Sp".encode())
            times.append(time.decode())

            # List converts bytes to list of bytes
            # The last two charachters are excluded - Newline char
            data_lines.append(list(data[:-2]))

        # The (-1/2) converts an unsigned byte to a real value
        power = np.single(np.asarray(data_lines) * (-1 / 2))
        times = np.double(np.asarray(times))

        return (power, times)


def tile_names():
    """List of reference and MWA tile names.

    Returns
    -------
    :return: list of tile names
    :rtype: list[str]

    """

    tiles = [
        "rf0XX",
        "rf0YY",
        "rf1XX",
        "rf1YY",
        "S06XX",
        "S06YY",
        "S07XX",
        "S07YY",
        "S08XX",
        "S08YY",
        "S09XX",
        "S09YY",
        "S10XX",
        "S10YY",
        "S12XX",
        "S12YY",
        "S29XX",
        "S29YY",
        "S30XX",
        "S30YY",
        "S31XX",
        "S31YY",
        "S32XX",
        "S32YY",
        "S33XX",
        "S33YY",
        "S34XX",
        "S34YY",
        "S35XX",
        "S35YY",
        "S36XX",
        "S36YY",
    ]

    return tiles


def tile_pairs(tiles):
    """
    Create a list of all possible AUT ref antenna pairs
    
    Parameters
    ----------
    :param tiles: list of tile names
    :type tiles: list[str]

    Returns
    -------
    :return: list of possible tile pairs
    :rtype: list[list, str]

    """

    tiles = tile_names()
    refs = tiles[:4]
    AUTS = tiles[4:]

    # Split the list into XX and YY lists
    refs_XX = refs[::2]
    refs_YY = refs[1::2]

    AUTS_XX = AUTS[::2]
    AUTS_YY = AUTS[1::2]

    # Create a list of pairs of tiles to be aligned.
    # Essentially all possible combinations of Refs, AUTS
    tile_pairs = []
    for pair in list(product(refs_XX, AUTS_XX)):
        tile_pairs.append(pair)

    for pair in list(product(refs_YY, AUTS_YY)):
        tile_pairs.append(pair)

    return tile_pairs


def time_tree(start_date, stop_date):
    """Split a date interval into 30 min observation chunks.

    This is used to travers the following directory tree.
    The data_root dir contains a sub-dir for every tile in
    :func:`~embers.condition_data.rf_data.tile_names` within
    which are dirs for every day, containg raw rf data files
    recorded every 30 minutes

    .. code-block:: text
    
        data_root
        ├── tile1
        │   └── YYYY-MM-DD
        │       └── tile1_YYYY-MM-DD-HH:MM.txt
        └── tile2 
            └── YYYY-MM-DD
                └── tile2_YYYY-MM-DD-HH:MM.txt

    Parameters
    ----------
    :param start_date: Start date in ``YYYY-MM-DD`` format
    :type start_date: str
    :param stop_date: Stop date is inclusive, in ``YYYY-MM-DD`` format. 
    :type stop_date: str

    Returns
    -------
    :returns:
        - dates - list of dates between `start_date` and `stop_date`
        - time_stamps - list of list of 30 minute observation in each day

    :rtype: list[str]


    """

    t_start = datetime.strptime(start_date, "%Y-%m-%d")
    t_stop = datetime.strptime(stop_date, "%Y-%m-%d")
    n_days = (t_stop - t_start).days

    dates = []
    time_stamps = []

    # Every Day
    for i in range(n_days + 1):
        day = t_start + timedelta(days=i)
        date = day.strftime("%Y-%m-%d")
        dates.append(date)
        d_t = []

        # Every 30 min in the day
        for j in range(48):
            t_delta = datetime.strptime(date, "%Y-%m-%d") + timedelta(minutes=30 * j)
            d_time = t_delta.strftime("%Y-%m-%d-%H:%M")
            d_t.append(d_time)

        time_stamps.append(d_t)

    return (dates, time_stamps)


def plt_waterfall(power, times, name):
    """
    Create waterfall `.pyplot.plot` object.

    waterfall created using parameters ``power``, ``times``
    from :func:`~embers.condition_data.rf_data.read_data`. 
    Default unix times are converted to a human readable
    *HH:MM* format.

    Parameters
    ----------
    :param power: `numpy.ndarray` object from :func:`~embers.condition_data.rf_data.read_data`
    :param times: `numpy.ndarray` object from :func:`~embers.condition_data.rf_data.read_data`
    :param name: name of rf data identifier to add to plot title

    Returns
    -------
        - plt - `.pyplot.plot` object

    """

    # setting dynamic range of waterfall to be 30 dB above the median
    power_median = np.median(power)
    image = power - power_median
    vmin = 0
    vmax = 30

    plt.style.use("dark_background")
    fig = plt.figure(figsize=(7, 10))
    ax = fig.add_axes([0.12, 0.1, 0.72, 0.85])
    im = ax.imshow(image, vmin=vmin, vmax=vmax, interpolation="none", cmap=_spec)
    cax = fig.add_axes([0.88, 0.1, 0.03, 0.85])
    fig.colorbar(im, cax=cax)
    ax.set_aspect("auto")
    ax.set_title(f"Waterfall plot: {name}")

    # Number of time steps on y-axis
    number_t = 5
    t_step = int(len(times) / (number_t - 1))
    times = list(times)
    times = times[::t_step]

    # Convert UNIX time to local HH:MM time
    t_tz = []
    for i in range(len(times)):
        perth_t = float(times[i]) + 28800  # 28800=+8GMT @ PERTH
        hms = time.strftime("%H:%M", time.gmtime(perth_t))
        t_tz.append(hms)

    # Frequency: x-axis
    start_freq = 137.15
    stop_freq = 138.55

    # X-axis stuff
    x_ax = image.shape[1]
    freqs = np.arange(start_freq, stop_freq, 0.25)
    x_ticks = np.arange(0, x_ax, (0.25 / 0.0125))  # .0125MHz/ch
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(freqs)
    ax.set_xlabel("Freq [MHz]")

    # Y-axis stuff
    y_ax = image.shape[0]
    y_ticks = np.arange(0, y_ax, t_step)
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(t_tz)
    ax.set_ylabel("Time [HH:MM]")

    return plt


def single_waterfall(rf_file, out_dir):
    """Save a waterfall plot from rf data file.
    
    Parameters
    ----------
    :param rf_file: path to a rf data file
    :type rf_file: str
    :param out_dir: path to output directory
    :type out_dir: str

    Returns
    -------
    :return: waterfall plot saved by `matplotlib.pyplot.savefig`

    """

    rf_name = Path(rf_file).stem

    power, times = read_data(rf_file)
    plt = plt_waterfall(power, times, rf_name)

    # Make out_dir if it doesn't exist
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    plt.savefig(f"{out_dir}/{rf_name}.png")
    plt.close()


def batch_waterfall(tile, time_stamp, data_dir, out_dir):
    """Save a waterfall plot for a batch of rf data files.

    Parameters
    ----------
    :param tile: tile name 
    :type tile: str
    :param time_stamp: when the rf observation began 
    :type time_stamp: str [*YYYY-MM-DD-HH:MM*]
    :param data_dir: path to root of data directory
    :type data_dir: str
    :param out_dir: path to output directory
    :type out_dir: str

    Returns
    -------
    :return: waterfall plot saved by `matplotlib.pyplot.savefig`
    
    Exeptions
    ---------
    :raises FileNotFoundError: input file does not exist

    """

    rf_name = f"{tile}_{time_stamp}"
    date = re.search(r"\d{4}.\d{2}.\d{2}", time_stamp)[0]
    rf_path = Path(f"{data_dir}/{tile}/{date}/{rf_name}.txt")

    try:
        open(rf_path, "r")
    except FileNotFoundError as e:
        return e
    except Exception as e:
        return e

    power, times = read_data(rf_path)
    plt = plt_waterfall(power, times, rf_name)

    # Make out_dir if it doesn't exist
    save_dir = Path(f"{out_dir}/waterfalls/{date}/{time_stamp}")
    save_dir.mkdir(parents=True, exist_ok=True)

    plt.savefig(f"{save_dir}/{rf_name}.png")
    plt.close()

    return f"Waterfall plot saved to {save_dir}/{rf_name}.png"
