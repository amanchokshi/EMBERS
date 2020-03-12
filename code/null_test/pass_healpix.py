import math
import json
import argparse
import numpy as np
import healpy as hp
import seaborn as sns
import concurrent.futures
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.stats import median_absolute_deviation as mad
from datetime import datetime, timedelta
from pathlib import Path

# Force matplotlib to not use X-Server backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
sys.path.append('../decode_rf_data')
sys.path.append('../sat_ephemeris')
sys.path.append('../sat_channels')
import rf_data as rf
from colormap import spectral
from sat_channels import time_tree, savgol_interp, time_filter
from channels_plt import plt_waterfall_pass, plt_channel, plt_hist, sat_plot

# Custom spectral colormap
cmap = spectral()
import sat_ids




def power_ephem(ref_file, chrono_file, sat_id, sat_chan):

    power, times = savgol_interp(ref_file, savgol_window, polyorder, interp_type, interp_freq )
    
    # To first order, let us consider the median to be the noise floor
    noise_f = np.median(power)
    noise_mad = mad(power, axis=None)
    
    # Scale the power to bring the median noise floor down to zero
    power = power - noise_f
    
    with open(chrono_file) as chrono:
        chrono_ephem = json.load(chrono)
    
        norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(len(chrono_ephem))]
        
        try:
            norad_index = norad_list.index(sat_id)
        
        except Exception:
            print(f'{sat_id}: not in {chrono_file}')
            return 0
        
        norad_ephem = chrono_ephem[norad_index]
            
        rise_ephem  = norad_ephem["time_array"][0] 
        set_ephem   = norad_ephem["time_array"][-1]
        
        intvl = time_filter(rise_ephem, set_ephem, np.asarray(times))
        
        if intvl != None:
            
            w_start, w_stop = intvl
        
            # Slice [crop] the power/times arrays to the times of sat pass
            channel_power = power[w_start:w_stop+1, sat_chan]
            times_c = times[w_start:w_stop+1]

    
            times_sat = norad_ephem["time_array"]
            
            # Crop ephem of all sats to size of norad_id sat
            intvl_ephem = time_filter(times_c[0], times_c[-1], np.asarray(times_sat))
            
            if intvl_ephem != None:
                e_0, e_1 = intvl_ephem
                    
                alt = norad_ephem["sat_alt"][e_0:e_1+1]
                az  = norad_ephem["sat_az"][e_0:e_1+1]

    return [channel_power, alt, az]

        
    

if __name__=='__main__':

    parser = argparse.ArgumentParser(description="""
        Project a sat pass onto a healpix map using ephemeris data
        """)
    
    parser.add_argument('--data_dir', metavar='\b', help='Dir where date is saved')
    #parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
    #parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/null_test/',help='Output directory. Default=./../../outputs/null_test/')
    parser.add_argument('--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
    parser.add_argument('--savgol_window', metavar='\b', default=151,help='Length of savgol window. Must be odd. Default=151')
    parser.add_argument('--polyorder', metavar='\b', default=1,help='Order of polynomial to fit to savgol window. Default=1')
    parser.add_argument('--interp_type', metavar='\b', default='cubic',help='Type of interpolation. Ex: cubic, linear, etc. Default=cubic')
    parser.add_argument('--interp_freq', metavar='\b', default=1,help='Frequency at which to resample smoothed data, in Hertz. Default=2')
    parser.add_argument('--nside', metavar='\b', default=32,help='Healpix Nside. Default = 32')
    
    args = parser.parse_args()
    
    data_dir =          args.data_dir
    chrono_dir =        args.chrono_dir
    #start_date =        args.start_date
    #stop_date =         args.stop_date
    out_dir =           args.out_dir
    savgol_window =     args.savgol_window
    polyorder =         args.polyorder
    interp_type =       args.interp_type
    interp_freq =       args.interp_freq
    nside =            args.nside
    
    # Save logs 
    #Path(out_dir).mkdir(parents=True, exist_ok=True)
    #sys.stdout = open(f'{out_dir}/Satellite_Channels_{start_date}_{stop_date}.txt', 'a')

    ref_file = Path(f'{data_dir}/rf0XX/2019-10-04/rf0XX_2019-10-04-16:00.txt')
    chrono_file = Path(f'{chrono_dir}/2019-10-04-16:00.json')

    sat_id = 41180
    sat_chan = 52
   
    sat_data = power_ephem(ref_file, chrono_file, sat_id, sat_chan)

    if sat_data != 0:
        channel_power, alt, az = sat_data
        print(len(channel_power), len(alt), len(az))



    ref_tile_map_dB_med  = [[] for pixel in range(hp.nside2npix(nside))]
    ref_tile_map_data_entries_counter=np.zeros(hp.nside2npix(nside))





