import math
import numpy as np

#TODO needed to add code dir to PYTHONPATH. Is this the best way?
from decode_rf_data.rf_data import read_data


from scipy import interpolate
from scipy.signal import savgol_filter


#tile_list = rf.tile_names()
def time_align(ref, tile):
    '''Time align power and timing arrays.
    
    This alignes the beginning and end of the time arrays by 
    minimising the difference between the start/stop times.
    It basically crops the array to only include
    overlapping times.

    Args:
        ref: path to reference data file
        tile: path to tile data file

    Returns:
        ref_t:      Reference time array
        ref_p:      Reference power array
        tile_t:     Tile time array
        tile_p:     Tile power array
    '''

    # Read time and power arrays from data files
    ref_p, ref_t = read_data(ref)
    tile_p, tile_t = read_data(tile)
    
    # Convert time arrays to numpy float arrays
    ref_t = np.asarray(ref_t).astype(float)
    tile_t = np.asarray(tile_t).astype(float)
    
    
    # Align the top of the time and power arrays
    delta_t = []
    if ref_t[0] >= tile_t[0]:
        for i in range(len(tile_t)):
            dt = (ref_t[0] - tile_t[i])
            delta_t.append(dt)
            if (dt <= 0 and abs(dt) >= delta_t[-1]):
                break
    
        abs_dt = np.absolute(delta_t)
        min_idx = (list(abs_dt).index(min(abs_dt)))
        
        tile_t = tile_t[min_idx:]
        tile_p = tile_p[min_idx:]
    
    else:
        for i in range(len(ref_t)):
            dt = (tile_t[0] - ref_t[i])
            delta_t.append(dt)
            if (dt <= 0 and abs(dt) >= delta_t[-1]):
                break
    
        abs_dt = np.absolute(delta_t)
        min_idx = (list(abs_dt).index(min(abs_dt)))
        
        ref_t = ref_t[min_idx:]
        ref_p = ref_p[min_idx:]
    
    
    delta_t = []
    if ref_t[-1] >= tile_t[-1]:
        for i in range(len(ref_t)):
            dt = (ref_t[-1 - i] - tile_t[-1])
            delta_t.append(dt)
            if (dt <= 0 and abs(dt) >= delta_t[-1]):
                break
    
        abs_dt = np.absolute(delta_t)
        min_idx = (list(abs_dt).index(min(abs_dt)))
        
        ref_t = ref_t[:-min_idx]
        ref_p = ref_p[:-min_idx]
    
    else:
        for i in range(len(tile_t)):
            dt = (tile_t[-1 - i] - ref_t[-1])
            delta_t.append(dt)
            if (dt <= 0 and abs(dt) >= delta_t[-1]):
                break
    
        abs_dt = np.absolute(delta_t)
        min_idx = (list(abs_dt).index(min(abs_dt)))
        
        tile_t = tile_t[:-min_idx]
        tile_p = tile_p[:-min_idx]
    
    return(ref_t, ref_p, tile_t, tile_p)
    
    

def savgol_interp(ref_t, ref_p, tile_t, tile_p, savgol_window =151, polyorder=1, interp_type='cubic', interp_freq=6):
    """Smooth and interpolate the power array

    Smooth the power arrays with a savgol filter
    and them interpolate to a given frequency.
    This makes the dimentions of the power arrays
    from thr reference antenna and the tile equal,
    allowing a one to one comparison at 
    corresponding times.

    
    Args:
        ref_t:          Reference time array
        ref_p:          Reference power array
        tile_t:         Tile time array
        tile_p:         Tile power array
        savgol_window:  Window size of savgol filer. Must be odd. Default = 151
        polyorder:      Order of polynomial to fit to savgol_window. Default = 1
        interp_type:    Type of interpolation. Ex: 'cubic', 'linear'. Default = cubic
        interp_freq:    The freqency to which power array is interpolated. Default = 6 Hz

    Returns:
        ref_p_aligned:  Aligned reference power array
        tile_p_aligned: Aligned tile power array
        time_array:     Time array corresponding to power arrays
    """


    print(len(tile_p[::, 1]))
    savgol_ref = savgol_filter(ref_p, savgol_window, polyorder, axis=0)
    savgol_tile = savgol_filter(tile_p, savgol_window, polyorder, axis=0)
    
    
    # Data is recorded at a range of frequencies, depending on the RF Explorer version
    # The Old models recoreded at 6-7Hz, while the new ones record at 8.5-9.2Hz.
    # To get around this, we interpolate the data to a desired freqeuncy.
    
    interp_freq = interp_freq #Hz
    
    # Using the start and stop times, we create an array of times at which to 
    # evaluate our interpolated data
    start_time = math.ceil(max(ref_t[0], tile_t[0]))
    stop_time = math.floor(min(ref_t[-1], tile_t[-1]))
    
    # Total length of observation in seconds
    time_seconds = stop_time - start_time
    
    # Array of times at which to evaluate the interpolated data
    time_array = np.linspace(start_time, stop_time, (time_seconds * interp_freq))
    
    # Interp1d output functions
    f = interpolate.interp1d(ref_t, savgol_ref, axis=0, kind=interp_type)
    g = interpolate.interp1d(tile_t, savgol_tile, axis=0, kind=interp_type)
    
    # New power array, evaluated at the desired frequency
    ref_p_aligned = f(time_array)
    tile_p_aligned = g(time_array)

    return (ref_p_aligned, tile_p_aligned, time_array)


if __name__ == '__main__':
    '''Plot a single channel to illustrate the savgol + interpolation
    
    Plots channel 59 of the sample date. This is the center 
    channel of the bright weather satellite seen in the bottom 
    right corner of the sample waterfall plots.'''
    
    import matplotlib
    # Enable x-window
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    
    ref_t, ref_p, tile_t, tile_p = time_align(
            './../../data/rf0XX_2019-10-10-02:30.txt',
            './../../data/S10XX_2019-10-10-02:30.txt'
            ) 

    ref_p_aligned, tile_p_aligned, time_array = savgol_interp(
            ref_t,
            ref_p,
            tile_t,
            tile_p
            )


    # Plots 
    plt.style.use('seaborn')

    plt.scatter(ref_t, ref_p[::, 59], color='#abcb89',marker='.', alpha=0.7, label='ref tile')
    plt.plot(time_array,ref_p_aligned[::, 59], color='#25a55f', alpha=0.9, label='ref savgol')	

    plt.scatter(tile_t, tile_p[::, 59], color='#ea7362',marker='.', alpha=0.7, label='aut tile')
    plt.plot(time_array,tile_p_aligned[::, 59], color='#a64942', alpha=0.9, label='aut savgol')

    leg = plt.legend(loc="upper left", frameon=True)
    leg.get_frame().set_facecolor('white')
    for l in leg.legendHandles:
        l.set_alpha(1)

    plt.xlabel('Time')
    plt.ylabel('Power')
    plt.tight_layout()
    plt.show()


