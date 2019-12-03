import math
import numpy as np

#TODO needed to add code dir to PYTHONPATH. Is this the best way?
from decode_rf_data.rf_data import read_data

import matplotlib
# Enable x-window
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from scipy import interpolate
from scipy.signal import savgol_filter


#tile_list = rf.tile_names()
def time_align(ref, tile):
    '''Time align power and timing arrays.
    
    This alignes the beginning and end of the time arrays by 
    minimising the difference between the start/stop times.
    It basically crops the array to only include
    overlapping times.
    '''

    ref_p, ref_t = read_data(ref)
    tile_p, tile_t = read_data(tile)
    
    ref_t = np.asarray(ref_t).astype(float)
    tile_t = np.asarray(tile_t).astype(float)
    
    
    delta_t = []
    if ref_t[0] >= tile_t[0]:
        for i in range(len(tile_t)):
            dt = (ref_t[0] - tile_t[i])
            delta_t.append(dt)
            if (dt <= 0 and abs(dt) >= delta_t[-1]):
                break
    
        abs_dt = np.absolute(delta_t)
        min_idx = (delta_t.index(min(abs_dt)))
        
        tile_t = tile_t[min_idx:]
        tile_p = tile_p[min_idx:]
    
    else:
        for i in range(len(ref_t)):
            dt = (tile_t[0] - ref_t[i])
            delta_t.append(dt)
            if (dt <= 0 and abs(dt) >= delta_t[-1]):
                break
    
        abs_dt = np.absolute(delta_t)
        min_idx = (delta_t.index(min(abs_dt)))
        
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
        min_idx = (delta_t.index(min(abs_dt)))
        
        ref_t = ref_t[:-min_idx]
        ref_p = ref_p[:-min_idx]
    
    else:
        for i in range(len(tile_t)):
            dt = (tile_t[-1 - i] - ref_t[-1])
            delta_t.append(dt)
            if (dt <= 0 and abs(dt) >= delta_t[-1]):
                break
    
        abs_dt = np.absolute(delta_t)
        min_idx = (delta_t.index(min(abs_dt)))
        
        tile_t = tile_t[:-min_idx]
        tile_p = tile_p[:-min_idx]
    
    return(ref_t, ref_p, tile_t, tile_p)
    
    
ref_t, ref_p, tile_t, tile_p = time_align('./../../data/rf0XX_2019-10-10-02:30.txt', './../../data/S10XX_2019-10-10-02:30.txt') 

#for channel in range(112):
for channel in list(range(56,61)):
    #savgol_sat = savgol_filter(ref_p, 123, 2, axis=0)
    
    savgol_ref = savgol_filter(ref_p[::, channel], 123, 1)
    savgol_tile = savgol_filter(tile_p[::, channel], 123, 1)
    
    
    # Data is recorded at a range of frequencies, depending on the RF Explorer version
    # The Old models recoreded at 6-7Hz, while the new ones record at 8.5-9.2Hz.
    # To get around this, we interpolate the data to a desired freqeuncy.
    
    interp_freq = 10 #Hz
    
    # Using the start and stop times, we create an array of times at which to 
    # evaluate our interpolated data
    start_time = math.ceil(max(ref_t[0], tile_t[0]))
    stop_time = math.floor(min(ref_t[-1], tile_t[-1]))
    
    # Total length of observation in seconds
    time_seconds = stop_time - start_time
    
    # Array of times at which to evaluate the interpolated data
    time_steps = np.linspace(start_time, stop_time, (time_seconds * interp_freq))
    
    # Interp1d output functions
    f = interpolate.interp1d(ref_t, savgol_ref, kind='cubic')
    g = interpolate.interp1d(tile_t, savgol_tile, kind='cubic')
    
    # New power array, evaluated at the desired frequency
    ref_p_n = f(time_steps)
    tile_p_n = g(time_steps)
    
    
    # Plots 
    plt.style.use('seaborn')
    
    plt.scatter(ref_t, ref_p[::, channel], color='#abcb89',marker='.', alpha=0.7, label='ref tile')
    plt.plot(ref_t,savgol_ref, color='#25a55f', alpha=0.9, label='ref savgol')
    
    
    plt.scatter(tile_t, tile_p[::, channel], color='#ea7362',marker='.', alpha=0.7, label='aut tile')
    plt.plot(tile_t,savgol_tile, color='#a64942', alpha=0.9, label='aut savgol')
    
    #plt.imshow(savgol_sat)
    #plt.axis('auto')
    
    leg = plt.legend(loc="upper left", frameon=True)
    leg.get_frame().set_facecolor('white')
    for l in leg.legendHandles:
        l.set_alpha(1)
    
    plt.xlabel('Time')
    plt.ylabel('Power')
    plt.xlim(1570646750, 1570647600)
    plt.title(f'Frequency Channel {channel}')
    plt.tight_layout()
    plt.savefig(f'channel-{channel}.png')
    #plt.show()
    plt.close()
