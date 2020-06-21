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
from embers.rf_tools.rf_data import read_data


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
    """Save an aligned set of rf data to an npz file. 

    A pair of rf data files are smoothed, interpolated and aligned
    with the :func:`~embers.condition_data.align_data.savgol_interp`.
    The output is written to a npz file and saved to an output 
    directory tree.

    Parameters
    ----------
    :param list[str] tile_pair: pair of ref and tile antenna names
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

    except FileNotFoundError as e:
        return e
    except Exception as e:
        return e
