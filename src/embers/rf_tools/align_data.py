"""
Align Data
----------
Tools to temporally align pairs of rf data files,
enabling comparisons between data sets

"""

import concurrent.futures
import logging
import math
import re
from itertools import repeat
from pathlib import Path

import numpy as np
from embers.rf_tools.rf_data import (read_data, tile_names, tile_pairs,
                                     time_tree)
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy.signal import savgol_filter


def savgol_interp(
    ref,
    tile,
    savgol_window_1=None,
    savgol_window_2=None,
    polyorder=None,
    interp_type=None,
    interp_freq=None,
):

    """Interpolate a power array followed by savgol smoothing.

    Interpolate to a given frequency,
    making the dimensions of the power arrays
    from reference and tile antennas equal,
    enabling comparisons between corresponding
    data points. Two level of savgol filter applied,
    first to capture deep nulls + small structure,
    and second level to smooth over noise.

    .. code-block:: python

        from embers.rf_tools.align_data import savgol_interp

         sg_interp_tuple = savgol_interp(
                            "~/embers-data/rf0XX.txt",
                            "~/embers-data/S06XX",
                            savgol_window_1=11,
                            savgol_window_2=15,
                            polyorder=2,
                            interp_type="cubic",
                            interp_freq=1)

        (ref_ali, tile_ali, time_array,
            ref_power, tile_power, ref_time, tile_time) = sg_interp_tuple

    :param ref: path to reference data file :class:`~str`
    :param tile: path to tile data file :class:`~str`
    :param savgol_window_1:  window size of savgol filer, must be odd :class:`~int`
    :param savgol_window_2:  window size of savgol filer, must be odd :class:`~int`
    :param polyorder: polynomial order to fit to savgol_window :class:`~int`
    :param interp_type: type of interpolation. Ex: 'cubic', 'linear' :class:`~str`
    :param interp_freq: freqency to which power array is interpolated in Hertz :class:`~int`

    :returns:
        A :class:`~tuple` (ref_ali, tile_ali, time_array, ref_power, tile_power, ref_time, tile_time)

        - ref_ali - aligned reference power array
        - tile_ali - aligned tile power array
        - time_array - corresponding to power arrays
        - ref_power - raw reference power array
        - tile_power - raw tile power array
        - ref_time - raw reference time array
        - tile_time - raw tile time array

    """

    # Read time and power arrays from data files
    ref_power, ref_time = read_data(ref)
    tile_power, tile_time = read_data(tile)

    # Round up/down to nearest integer of time
    start_time = math.ceil(max(ref_time[0], tile_time[0]))
    stop_time = math.floor(min(ref_time[-1], tile_time[-1]))

    # Array of times at which to evaluate the interpolated data
    time_array = np.arange(start_time, stop_time, (1 / interp_freq))

    # Mathematical interpolation functions
    f = interpolate.interp1d(ref_time, ref_power, axis=0, kind=interp_type)
    g = interpolate.interp1d(tile_time, tile_power, axis=0, kind=interp_type)

    # New power array, evaluated at the desired frequency
    ref_ali = f(time_array)
    tile_ali = g(time_array)

    # Savgol level 1. Capture nulls / small scale structure
    ref_ali = savgol_filter(ref_ali, savgol_window_1, polyorder, axis=0)
    tile_ali = savgol_filter(tile_ali, savgol_window_1, polyorder, axis=0)

    # Savgol level 2. Smooth noise
    ref_ali = savgol_filter(ref_ali, savgol_window_2, polyorder, axis=0)
    tile_ali = savgol_filter(tile_ali, savgol_window_2, polyorder, axis=0)

    return (ref_ali, tile_ali, time_array, ref_power, tile_power, ref_time, tile_time)


def plot_savgol_interp(
    ref=None,
    tile=None,
    savgol_window_1=None,
    savgol_window_2=None,
    polyorder=None,
    interp_type=None,
    interp_freq=None,
    channel=None,
    out_dir=None,
):
    """Plot single channel of power arrays to visualise :func:`~embers.rf_tools.align_data.savgol_interp`.

    Create a plot of a single channel of raw :samp:`rf_data` from reference and tile power arrays, along
    with the outputs of :func:`~embers.rf_tools.align_data.savgol_interp` to visualise the effects of
    interpolation and savgol smoothing.

    :param ref: path to reference data file :class:`~str`
    :param tile: path to tile data file :class:`~str`
    :param savgol_window_1:  window size of savgol filer, must be odd :class:`~int`
    :param savgol_window_2:  window size of savgol filer, must be odd :class:`~int`
    :param polyorder: polynomial order to fit to savgol_window :class:`~int`
    :param interp_type: type of interpolation. Ex: 'cubic', 'linear' :class:`~str`
    :param interp_freq: freqency to which power array is interpolated in Hertz :class:`~int`
    :param channel: index of single frequency channel :class:`~int`
    :param out_dir: path to output directory :class:`~str`

    :returns:
        single freqency savgol_interp plot saved to :samp:`out_dir`

    """

    (
        ref_ali,
        tile_ali,
        time_array,
        ref_power,
        tile_power,
        ref_time,
        tile_time,
    ) = savgol_interp(
        ref=ref,
        tile=tile,
        savgol_window_1=savgol_window_1,
        savgol_window_2=savgol_window_2,
        polyorder=polyorder,
        interp_type=interp_type,
        interp_freq=interp_freq,
    )

    # Sample align plot
    plt.style.use("seaborn")
    plt.rcParams["figure.figsize"] = (9, 6)

    # convert times to minuts from first datapoint
    time_array = (time_array - time_array[0]) / 60
    ref_time = (ref_time - ref_time[0]) / 60
    tile_time = (tile_time - tile_time[0]) / 60

    plt.plot(
        time_array,
        tile_ali[::, channel],
        color="#e23a4e",
        alpha=0.9,
        label="tile savgol",
    )
    plt.scatter(
        tile_time,
        tile_power[::, channel],
        color="#f78b51",
        marker=".",
        alpha=0.6,
        label="tile raw",
    )
    plt.plot(
        time_array,
        ref_ali[::, channel],
        color="#252b40",
        alpha=0.9,
        label="ref savgol",
    )
    plt.scatter(
        ref_time,
        ref_power[::, channel],
        color="#6a82bb",
        marker=".",
        alpha=0.6,
        label="ref raw",
    )

    leg = plt.legend(loc="upper left", frameon=True)
    leg.get_frame().set_facecolor("white")
    for leg in leg.legendHandles:
        leg.set_alpha(1)

    plt.ylim(-110, -20)
    plt.ylabel("Raw Power [dBm]")
    plt.xlabel("Time [min]")
    plt.tight_layout()
    Path(f"{out_dir}").mkdir(parents=True, exist_ok=True)
    plt.savefig(f"{out_dir}/savgol_interp_sample.png")


def save_aligned(
    tile_pair,
    time_stamp,
    savgol_window_1,
    savgol_window_2,
    polyorder,
    interp_type,
    interp_freq,
    data_dir,
    out_dir,
):
    """Save an aligned set of rf data with :func:`~numpy.savez_compressed` to an :samp:`npz` file.

    A pair of rf data files are smoothed, interpolated and aligned
    with the :func:`~embers.rf_tools.align_data.savgol_interp`.
    with the output written to a :samp:`npz` file and saved to an output
    directory tree.

    .. code-block:: python

        from embers.rf_tools.align_data import save_aligned

        savgol_interp(
                    ["rf0XX", "S06XX"],
                    "2020-01-01-00:00"
                    savgol_window_1=11,
                    savgol_window_2=15,
                    polyorder=2,
                    interp_type="cubic",
                    interp_freq=1,
                    "~/embers-data/",
                    "~/embers-outputs")


    :param tile_pair: pair of ref and tile antenna names from :func:`~embers.rf_tools.rf_data.tile_pairs` :class:`list`
    :param time_stamp: time when rf observation began. In YYYY-MM-DD-HH-MM format :class:`~str`
    :param savgol_window_1:  window size of savgol filer, must be odd :class:`~int`
    :param savgol_window_2:  window size of savgol filer, must be odd :class:`~int`
    :param polyorder: polynomial order to fit to savgol_window :class:`~int`
    :param interp_type: type of interpolation. Ex: 'cubic', 'linear' :class:`~str`
    :param interp_freq: freqency to which power array is interpolated :class:`~int`
    :param data_dir: root of data dir where rf data is located :class:`~str`
    :param out_dir: relative path to output directory :class:`~str`

    :return:
        - aligned rf data saved to :samp:`npz` file by :func:`~numpy.savez_compressed`

    :raises FileNotFoundError: an input file does not exist

    """

    date = re.search(r"\d{4}.\d{2}.\d{2}", time_stamp)[0]

    ref = tile_pair[0]
    tile = tile_pair[1]
    ref_file = f"{data_dir}/{ref}/{date}/{ref}_{time_stamp}.txt"
    tile_file = f"{data_dir}/{tile}/{date}/{tile}_{time_stamp}.txt"

    try:
        ref_ali, tile_ali, time_array, _, _, _, _ = savgol_interp(
            ref_file,
            tile_file,
            savgol_window_1=savgol_window_1,
            savgol_window_2=savgol_window_2,
            polyorder=polyorder,
            interp_type=interp_type,
            interp_freq=interp_freq,
        )

        # creates output directory if it doesn't exist
        save_dir = Path(f"{out_dir}/{date}/{time_stamp}")
        save_dir.mkdir(parents=True, exist_ok=True)

        # Convert the power array to float32
        # Convert list of times to float64 (double)
        # Save as compressed npz file. Seems to drastically reduce size
        np.savez_compressed(
            f"{save_dir}/{ref}_{tile}_{time_stamp}_aligned.npz",
            ref_ali=np.single(ref_ali),
            tile_ali=np.single(tile_ali),
            time_array=np.double(time_array),
        )

        return f"Saved aligned file to {save_dir}/{ref}_{tile}_{time_stamp}_aligned.npz"

    except Exception as e:
        return e


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
    max_cores=None,
):
    """Temporally align all RF files within a date interval using :func:`~embers.rf_tools.align_data.save_aligned`.


    :param start_date: In YYYY-MM-DD format :class:`~str`
    :param stop_date: In YYYY-MM-DD format :class:`~str`
    :param savgol_window_1:  window size of savgol filer, must be odd :class:`~int`
    :param savgol_window_2:  window size of savgol filer, must be odd :class:`~int`
    :param polyorder: polynomial order to fit to savgol_window :class:`~int`
    :param interp_type: type of interpolation. Ex: 'cubic', 'linear' :class:`~str`
    :param interp_freq: freqency to which power array is interpolated :class:`~int`
    :param data_dir: root of data dir where rf data is located :class:`~str`
    :param out_dir: relative path to output directory :class:`~str`
    :param max_cores: Maximum number of cores to be used by this script. Default=None, which means that all available cores are used

    :return:
        - aligned rf data saved to :samp:`npz` file by :func:`~numpy.savez_compressed` in :samp:`out_dir`

    """

    dates, time_stamps = time_tree(start_date, stop_date)

    # Logging config
    log_dir = Path(f"{out_dir}")
    log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=f"{out_dir}/align_batch.log",
        level=logging.INFO,
        format="%(levelname)s: %(funcName)s: %(message)s",
    )

    for pair in tile_pairs(tile_names()):

        for day in range(len(dates)):

            with concurrent.futures.ProcessPoolExecutor(
                max_workers=max_cores
            ) as executor:
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
