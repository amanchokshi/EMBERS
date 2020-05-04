import math
import time
import numpy as np

import sys
sys.path.append('../decode_rf_data')
from rf_data import read_data

from scipy import interpolate
from scipy.signal import savgol_filter


def savgol_interp(ref, tile, savgol_window_1 =None,savgol_window_2=None, polyorder=None, interp_type=None, interp_freq=None):
    """Smooth and interpolate the power array

    Smooth the power arrays with a savgol filter
    and them interpolate to a given frequency.
    This makes the dimentions of the power arrays
    from thr reference antenna and the tile equal,
    allowing a one to one comparison at 
    corresponding times.

    
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
    
    # Read time and power arrays from data files
    ref_p, ref_t = read_data(ref)
    tile_p, tile_t = read_data(tile)

    #savgol_ref = savgol_filter(ref_p, savgol_window, polyorder, axis=0)
    #savgol_tile = savgol_filter(tile_p, savgol_window, polyorder, axis=0)

    
    # Data is recorded at a range of frequencies, depending on the RF Explorer version
    # The Old models recoreded at 6-7Hz, while the new ones record at 8.5-9.2Hz.
    # To get around this, we interpolate the data to a desired freqeuncy.
    
    interp_freq = interp_freq #Hz
    
    # Using the start and stop times, we create an array of times at which to 
    # evaluate our interpolated data
    # Round up/down to nearest integer of time
    start_time = math.ceil(max(ref_t[0], tile_t[0]))
    stop_time = math.floor(min(ref_t[-1], tile_t[-1]))
    
    # Array of times at which to evaluate the interpolated data
    time_array = np.arange(start_time, stop_time, (1 / interp_freq))
    
    # Interp1d output functions
    #f = interpolate.interp1d(ref_t, savgol_ref, axis=0, kind=interp_type)
    #g = interpolate.interp1d(tile_t, savgol_tile, axis=0, kind=interp_type)
    
    f = interpolate.interp1d(ref_t, ref_p, axis=0, kind=interp_type)
    g = interpolate.interp1d(tile_t, tile_p, axis=0, kind=interp_type)
    
    # New power array, evaluated at the desired frequency
    ref_p_aligned = f(time_array)
    tile_p_aligned = g(time_array)
    
    ref_p_aligned = savgol_filter(ref_p_aligned, savgol_window_1, polyorder, axis=0)
    tile_p_aligned = savgol_filter(tile_p_aligned, savgol_window_1, polyorder, axis=0)
    
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
            './../../data/rf0XX_2019-10-10-02:30.txt',
            './../../data/S10XX_2019-10-10-02:30.txt',
            savgol_window_1 =11,
            savgol_window_2 =15,
            polyorder=2,
            interp_type='cubic',
            interp_freq=1
            )


#    ch =23
#    ref_t, ref_p, tile_t, tile_p, ref_p_aligned, tile_p_aligned, time_array = savgol_interp(
#            '../../../tiles_data/rf0XX/2019-09-12/rf0XX_2019-09-12-09:30.txt',
#            '../../../tiles_data/S06XX/2019-09-12/S06XX_2019-09-12-09:30.txt',
#            savgol_window_1 =11,
#            savgol_window_2 =15,
#            polyorder=2,
#            interp_type='cubic',
#            interp_freq=1
#            )

#    ch =13
#    ref_t, ref_p, tile_t, tile_p, ref_p_aligned, tile_p_aligned, time_array = savgol_interp(
#            '../../../tiles_data/rf0XX/2019-09-14/rf0XX_2019-09-14-11:30.txt',
#            '../../../tiles_data/S06XX/2019-09-14/S06XX_2019-09-14-11:30.txt',
#            savgol_window_1 =11,
#            savgol_window_2 =15,
#            polyorder=2,
#            interp_type='cubic',
#            interp_freq=1
#            )


#    ch =8
#    ref_t, ref_p, tile_t, tile_p, ref_p_aligned, tile_p_aligned, time_array = savgol_interp(
#            '../../../tiles_data/rf0XX/2019-09-15/rf0XX_2019-09-15-11:00.txt',
#            '../../../tiles_data/S06XX/2019-09-15/S06XX_2019-09-15-11:00.txt',
#            savgol_window_1 =11,
#            savgol_window_2 =15,
#            polyorder=2,
#            interp_type='cubic',
#            interp_freq=1
#            )

# Plots
    fig = plt.figure(figsize=(12,9))
    plt.style.use('seaborn')
#    plt.rcParams["figure.figsize"] = (6,9)
#    plt.rcParams.update({
#    "lines.color": "white",
#    "text.color": "black",
#    "axes.facecolor": "#b3b3b3",
#    "axes.edgecolor": "lightgray",
#    "axes.labelcolor": "white",
#    "xtick.color": "white",
#    "ytick.color": "white",
#    "grid.color": "lightgray",
#    "figure.facecolor": "#2F1543",
#    "figure.edgecolor": "black",
#    "savefig.facecolor": "#2F1543",
#    "savefig.edgecolor": "black"})


    plt.plot(time_array,tile_p_aligned[::, ch], color='#e23a4e', alpha=0.9, label='aut savgol')
    #plt.fill_between(time_array,tile_p_aligned[::, 59], np.full(len(time_array), min(ref_p_aligned[::, 59])), color='#ff5453', alpha=0.8)	
    plt.scatter(tile_t, tile_p[::, ch], color='#f78b51',marker='.', alpha=0.6, label='aut tile')
    
    plt.plot(time_array,ref_p_aligned[::, ch], color='#252b40', alpha=0.9, label='ref savgol')	
    #plt.fill_between(time_array,ref_p_aligned[::, 59], np.full(len(time_array), min(ref_p_aligned[::, 59])), color='#252b40', alpha=0.8)	
    plt.scatter(ref_t, ref_p[::, ch], color='#6a82bb',marker='.', alpha=0.6, label='ref tile')
#    plt.scatter(ref_t, ref_p[::, 59], color='#fe6845',marker='.', alpha=0.7)
#    plt.plot(time_array,ref_p_aligned[::, 59], color='#ec9b3b', alpha=0.9, label='Ref Tile')	
#
#    plt.scatter(tile_t, tile_p[::, 59], color='#3a1f5d',marker='.', alpha=0.7)
#    plt.plot(time_array,tile_p_aligned[::, 59], color='#916dd5', alpha=0.9, label='MWA Tile')

    leg = plt.legend(loc="upper left", frameon=True)
    leg.get_frame().set_facecolor('white')
    for l in leg.legendHandles:
        l.set_alpha(1)

   

#    t_sta = 1570646700
#    t_sto = 1570647550
#    plt.xlim(t_sta, t_sto)
#    
#    times = np.linspace(t_sta, t_sto, 5)
#    
#    t_tz = []
#    
#    # Convert UNIX time to local HH:MM time
#    for i in range(len(times)):
#        
#        perth_t = float(times[i])+28800 #28800=+8GMT @ PERTH
#        hms = time.strftime('%H:%M', time.gmtime(perth_t))
#        t_tz.append(hms)


    
    #plt.xlabel('Time [HH:MM]')
    #plt.xticks(times, t_tz)
    plt.ylim(-120, 5)
    plt.ylabel('Power')
    plt.tight_layout()
    plt.savefig('test.png', dpi=300)
    #plt.show()


