import sys
sys.path.append('../decode_rf_data')

import numpy as np
from rf_data import read_data

import matplotlib
# Enable x-window display
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


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

print(f'△ T = {(ref_t[-1] - ref_t[0]) - (tile_t[-1] - tile_t[0])} seconds')


plt.style.use('seaborn')
plt.plot(ref_p[::, 4], color='#7da87b', alpha=0.9)
#plt.plot(tile_p[:, 4], color='#ed6663', alpha=0.8)
plt.tight_layout()
plt.show()
