"""
Align Data
==========
Tools to time align pairs of rf data files,
enabling comparisons between corresponding 
data points.

"""

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
        - ref_power - raw reference power array
        - tile_power - raw tile power array
        - time_array - time array corresponding to power arrays

    :rtype: (float, numpy.array(float))

    """

    # Read time and power arrays from data files
    ref_power, ref_t = read_data(ref)
    tile_power, tile_t = read_data(tile)

    # Round up/down to nearest integer of time
    start_time = math.ceil(max(ref_t[0], tile_t[0]))
    stop_time = math.floor(min(ref_t[-1], tile_t[-1]))

    # Array of times at which to evaluate the interpolated data
    time_array = np.arange(start_time, stop_time, (1 / interp_freq))

    # Mathematical interpolation functions
    f = interpolate.interp1d(ref_t, ref_power, axis=0, kind=interp_type)
    g = interpolate.interp1d(tile_t, tile_power, axis=0, kind=interp_type)

    # New power array, evaluated at the desired frequency
    ref_ali = f(time_array)
    tile_ali = g(time_array)

    # Savgol level 1. Capture nulls / small scale structure
    ref_ali = savgol_filter(ref_ali, savgol_window_1, polyorder, axis=0)
    tile_ali = savgol_filter(tile_ali, savgol_window_1, polyorder, axis=0)

    # Savgol level 2. Smooth noise
    ref_ali = savgol_filter(ref_ali, savgol_window_2, polyorder, axis=0)
    tile_ali = savgol_filter(tile_ali, savgol_window_2, polyorder, axis=0)

    return (ref_ali, tile_ali, ref_power, tile_power, time_array)


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
