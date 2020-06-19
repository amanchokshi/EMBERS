"""
Align Data
==========
Tools to time align pairs of rf data files,
enabling comparisons between corresponding 
data points.

"""

import re
import math
import time
import numpy as np
from pathlib import Path
from scipy import interpolate
from scipy.signal import savgol_filter
from embers.condition_data.rf_data import read_data


def savgol_interp(
    ref,
    tile,
    savgol_window_1=None,
    savgol_window_2=None,
    polyorder=None,
    interp_type=None,
    interp_freq=None,
):
    """Interpolate the power array followed by savgol smoothing.

    Interpolate to a given frequency which 
    makes the dimensions of the power arrays
    from the reference and tile antennas equal,
    enabling comparisons between corresponding
    data points. Two level of savgol filter applied, 
    first to capture deep nulls + small structure, 
    and second level to smooth over noise.

    
    Parameters
    ----------
    :param str ref: path to reference data file
    :param str tile: path to tile data file
    :param int savgol_window_1:  window size of savgol filer, must be odd
    :param int savgol_window_2:  window size of savgol filer, must be odd
    :param int polyorder: polynomial order to fit to savgol_window
    :param str interp_type: type of interpolation. Ex: 'cubic', 'linear'
    :param int interp_freq: freqency to which power array is interpolated

    Returns
    -------
    :returns:
        - ref_ali - aligned reference power array
        - tile_ali - aligned tile power array
        - time_array - time array corresponding to power arrays
        - ref_power - raw reference power array
        - tile_power - raw tile power array
        - ref_time - raw reference time array
        - tile_time - raw tile time array

    :rtype: (float, numpy.array(float))

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


def save_aligned(
    ref,
    aut,
    time_stamp,
    savgol_window_1,
    savgol_window_2,
    polyorder,
    interp_type,
    interp_freq,
    data_dir,
    out_dir,
):
    """Save an aligned set of rf data to an npz file. 

    A pair of rf data files are smoothed, interpolated and aligned
    with the :func:`~embers.condition_data.align_data.savgol_interp`.
    The output is written to a npz file and saved to an output 
    directory tree.

    Parameters
    ----------
    :param str ref: name of reference antenna
    :param str tile: name of tile antenna
    :param str time_stamp: time when rf observation began. In YYYY-MM-DD-HH-MM format
    :param int savgol_window_1:  window size of savgol filer, must be odd
    :param int savgol_window_2:  window size of savgol filer, must be odd
    :param int polyorder: polynomial order to fit to savgol_window
    :param str interp_type: type of interpolation. Ex: 'cubic', 'linear'
    :param int interp_freq: freqency to which power array is interpolated
    :param str data_dir: root of data dir where rf data is located
    :param str out_dir: relative path to output directory

    Returns
    -------
    :return: aligned rf data saved by `numpy.savez_compressed`
    
    Exeptions
    ---------
    :raises FileNotFoundError: an input file does not exist

    """

    date = re.search(r"\d{4}.\d{2}.\d{2}", time_stamp)[0]

    ref_file = f"{data_dir}/{ref}/{date}/{ref}_{time_stamp}.txt"
    tile_file = f"{data_dir}/{tile}/{date}/{tile}_{time_stamp}.txt"

    if Path(ref_file) and Path(tile_file):

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
            f"{save_dir}/{ref}_{aut}_{time_stamp}_aligned.npz",
            ref_ali=np.single(ref_ali),
            tile_ali=np.single(tile_ali),
            time_array=np.double(time_array),
        )

        return f"Saved aligned file to {save_dir}/{ref}_{aut}_{time_stamp}_aligned.npz"

    else:
        return f"FileNotFoundError: either {ref_file} or {tile_file} missing"


if __name__ == "__main__":
    """Plot a single channel to illustrate the savgol + interpolation
    
    Plots channel 59 of the sample date. This is the center 
    channel of the bright weather satellite seen in the bottom 
    right corner of the sample waterfall plots."""

    import matplotlib

    # Enable x-window
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ch = 59
    (
        ref_t,
        ref_p,
        tile_t,
        tile_p,
        ref_p_aligned,
        tile_p_aligned,
        time_array,
    ) = savgol_interp(
        "./../../data/rf_data/rf0XX/2019-10-10/rf0XX_2019-10-10-02:30.txt",
        "./../../data/rf_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt",
        savgol_window_1=11,
        savgol_window_2=15,
        polyorder=2,
        interp_type="cubic",
        interp_freq=1,
    )

    # Sample align plot
    plt.style.use("seaborn")
    plt.rcParams["figure.figsize"] = (9, 6)

    time_array = (time_array - time_array[0]) / 60
    ref_t = (ref_t - ref_t[0]) / 60
    tile_t = (tile_t - tile_t[0]) / 60

    plt.plot(
        time_array,
        tile_p_aligned[::, ch],
        color="#e23a4e",
        alpha=0.9,
        label="aut savgol",
    )
    plt.scatter(
        tile_t, tile_p[::, ch], color="#f78b51", marker=".", alpha=0.6, label="aut tile"
    )

    plt.plot(
        time_array,
        ref_p_aligned[::, ch],
        color="#252b40",
        alpha=0.9,
        label="ref savgol",
    )
    plt.scatter(
        ref_t, ref_p[::, ch], color="#6a82bb", marker=".", alpha=0.6, label="ref tile"
    )

    leg = plt.legend(loc="upper left", frameon=True)
    leg.get_frame().set_facecolor("white")
    for l in leg.legendHandles:
        l.set_alpha(1)

    plt.ylim(-110, -20)
    plt.ylabel("Power [dBm]")
    plt.xlabel("Time [min]")
    plt.tight_layout()
    Path("../../outputs/align_data/").mkdir(parents=True, exist_ok=True)
    plt.savefig("../../outputs/align_data/align_data_test.png", dpi=300)
