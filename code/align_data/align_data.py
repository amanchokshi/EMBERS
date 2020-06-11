import math
import time
import numpy as np
from pathlib import Path
from scipy import interpolate
from scipy.signal import savgol_filter


def savgol_interp(ref, tile, savgol_window_1 =None,savgol_window_2=None, polyorder=None, interp_type=None, interp_freq=None):
    """Smooth and interpolate the power array

    Interpolate to a given frequency which 
    makes the dimentions of the power arrays
    from the reference antenna and the tile equal,
    allowing a one to one comparison at corresponding times.

    Two level of savgol filter applied, first to capture
    deep nulls + small structure, and second level to
    smooth over noise.

    
    Args:
        ref:            Path to reference data file
        tile:           Path to tile data file
        savgol_window_1:  Window size of savgol filer. Must be odd. Default = 151
        savgol_window_2:  Window size of savgol filer. Must be odd. Default = 151
        polyorder:      Order of polynomial to fit to savgol_window. Default = 1
        interp_type:    Type of interpolation. Ex: 'cubic', 'linear'. Default = cubic
        interp_freq:    The freqency to which power array is interpolated. Default = 6 Hz

    Returns:
        ref_p_aligned:  Aligned reference power array
        tile_p_aligned: Aligned tile power array
        time_array:     Time array corresponding to power arrays
    """

    # Import custom rf_data module
    import sys
    sys.path.append('../decode_rf_data')
    from rf_data import read_data
    
    # Read time and power arrays from data files
    ref_p, ref_t = read_data(ref)
    tile_p, tile_t = read_data(tile)

    # Data is recorded at a range of frequencies, depending on the RF Explorer version
    # The Old models recoreded at 6-7Hz, while the new ones record at 8.5-9.2Hz.
    # To get around this, we interpolate the data to a desired freqeuncy.
    
    
    # Using the start and stop times, we create an array of times at which to 
    # evaluate our interpolated data
    # Round up/down to nearest integer of time
    start_time = math.ceil(max(ref_t[0], tile_t[0]))
    stop_time = math.floor(min(ref_t[-1], tile_t[-1]))
    
    # Array of times at which to evaluate the interpolated data
    time_array = np.arange(start_time, stop_time, (1 / interp_freq))
    
    f = interpolate.interp1d(ref_t, ref_p, axis=0, kind=interp_type)
    g = interpolate.interp1d(tile_t, tile_p, axis=0, kind=interp_type)
    
    # New power array, evaluated at the desired frequency
    ref_p_aligned = f(time_array)
    tile_p_aligned = g(time_array)
   
    # Savgol level 1. Capture nulls / small scale structure
    ref_p_aligned = savgol_filter(ref_p_aligned, savgol_window_1, polyorder, axis=0)
    tile_p_aligned = savgol_filter(tile_p_aligned, savgol_window_1, polyorder, axis=0)
    
    # Savgol level 2. Smooth noise
    ref_p_aligned = savgol_filter(ref_p_aligned,   savgol_window_2, polyorder, axis=0)
    tile_p_aligned = savgol_filter(tile_p_aligned, savgol_window_2, polyorder, axis=0)
    
    return (ref_t, ref_p, tile_t, tile_p, ref_p_aligned, tile_p_aligned, time_array)


if __name__ == '__main__':
    '''Plot a single channel to illustrate the savgol + interpolation
    
    Plots channel 59 of the sample date. This is the center 
    channel of the bright weather satellite seen in the bottom 
    right corner of the sample waterfall plots.'''
    
    import matplotlib
    # Enable x-window
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    

    ch =59
    ref_t, ref_p, tile_t, tile_p, ref_p_aligned, tile_p_aligned, time_array = savgol_interp(
            './../../data/rf_data/rf0XX/2019-10-10/rf0XX_2019-10-10-02:30.txt',
            './../../data/rf_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt',
            savgol_window_1 =11,
            savgol_window_2 =15,
            polyorder=2,
            interp_type='cubic',
            interp_freq=1
            )



    # Sample align plot
    plt.style.use('seaborn')
    plt.rcParams["figure.figsize"] = (9,6)
    
    time_array = (time_array - time_array[0])/60
    ref_t = (ref_t - ref_t[0])/60
    tile_t = (tile_t - tile_t[0])/60


    plt.plot(time_array,tile_p_aligned[::, ch], color='#e23a4e', alpha=0.9, label='aut savgol')
    plt.scatter(tile_t, tile_p[::, ch], color='#f78b51',marker='.', alpha=0.6, label='aut tile')
    
    plt.plot(time_array,ref_p_aligned[::, ch], color='#252b40', alpha=0.9, label='ref savgol')	
    plt.scatter(ref_t, ref_p[::, ch], color='#6a82bb',marker='.', alpha=0.6, label='ref tile')

    leg = plt.legend(loc="upper left", frameon=True)
    leg.get_frame().set_facecolor('white')
    for l in leg.legendHandles:
        l.set_alpha(1)

   
    plt.ylim(-110, -20)
    plt.ylabel('Power [dBm]')
    plt.xlabel('Time [min]')
    plt.tight_layout()
    Path('../../outputs/align_data/').mkdir(parents=True, exist_ok=True)
    plt.savefig('../../outputs/align_data/align_data_test.png', dpi=300)


