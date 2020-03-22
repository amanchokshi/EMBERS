import math
import json
import numpy as np
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.stats import median_absolute_deviation as mad
from datetime import datetime, timedelta
from pathlib import Path

import sys
sys.path.append('../decode_rf_data')
sys.path.append('../sat_ephemeris')
import rf_data as rf
import sat_ids

## Custom spectral colormap
from colormap import spectral
cmap = spectral()


def time_tree(start_date, stop_date):
    '''Split the time interval into 30 min obs'''

    t_start = datetime.strptime(start_date, '%Y-%m-%d')
    t_stop = datetime.strptime(stop_date, '%Y-%m-%d')
    n_days = (t_stop - t_start).days
    
    dates = []
    date_time = []
   
    # Every Day
    for i in range(n_days+1):
        day = t_start + timedelta(days=i)
        date = day.strftime('%Y-%m-%d')
        dates.append(date)
        d_t = []
        
        # Every 30 min in the day
        for j in range(48):
            t_delta = datetime.strptime(date,'%Y-%m-%d') + timedelta(minutes=30*j)
            d_time = t_delta.strftime('%Y-%m-%d-%H:%M')
            d_t.append(d_time)
    
        date_time.append(d_t)

    return (dates, date_time)


def savgol_interp(ref, savgol_window =None, polyorder=None, interp_type=None, interp_freq=None):
    """Smooth and interpolate the power array with a savgol filter
    
    Args:
        ref:            Path to reference data file
        savgol_window:  Window size of savgol filer. Must be odd. Default = 151
        polyorder:      Order of polynomial to fit to savgol_window. Default = 1
        interp_type:    Type of interpolation. Ex: 'cubic', 'linear'. Default = cubic
        interp_freq:    The freqency to which power array is interpolated. Default = 6 Hz

    Returns:
        power_smooth:   Aligned reference power array
        time_smooth:    Time array corresponding to power arrays
    """
    
    
    power, times = rf.read_data(ref)
    savgol_pow = savgol_filter(power, savgol_window, polyorder, axis=0)
    time_smooth = np.arange(math.ceil(times[0]), math.floor(times[-1]), (1 / interp_freq))
    ref_sav_interp = interpolate.interp1d(times, savgol_pow, axis=0, kind=interp_type)
    power_smooth = ref_sav_interp(time_smooth)
    
    return (power_smooth, time_smooth)


def center_of_gravity(channel_power, times_c):
    '''Determine center of gravity of channel power
    
    Args:
        channel_power:  Power in one channel
        times_c:        Corresponding time array
    '''

    # a list of indices for a column, and the center of this list
    #t = np.arange(channel_power.shape[0])
    center = times_c[(len(times_c) - 1)//2]
   
    # We determine the center of gravity for each of the possible sat channels
    cog = np.round(np.sum(channel_power * times_c) / np.sum(channel_power), 1)
    
    # Delta c. distance b/w cog and center
    del_c = np.absolute(cog - center)

    # Fractional offset from center
    frac_cen_offset = del_c/(channel_power.shape[0])

    return (center, cog, frac_cen_offset)


def time_filter(s_rise, s_set, times):
    '''Slice obs window to size of ephem of norad sat
    Args:
        s_rise          : ephem rise time
        s_set           : ephem set time
        times           : time array of obs
        '''
    # I. sat rises before times, sets within times window
    if (s_rise < times[0] and s_set > times[0] and s_set <= times[-1]):
        i_0 = np.where(times == times[0])[0][0]
        i_1 = np.where(times == s_set)[0][0]
        intvl = [i_0, i_1]
    
    # II. sat rises and sets within times
    elif (s_rise >= times[0] and s_set <= times[-1]):
        i_0 = np.where(times == s_rise)[0][0]
        i_1 = np.where(times == s_set)[0][0]
        intvl = [i_0, i_1]
    
    # III. sat rises within times, and sets after
    elif (s_rise >= times[0] and s_rise < times[-1] and s_set > times[-1]):
        i_0 = np.where(times == s_rise)[0][0]
        i_1 = np.where(times == times[-1])[0][0]
        intvl = [i_0, i_1]
    
    # IV. sat rises before times and sets after
    elif (s_rise < times[0] and s_set > times[-1]):
        i_0 = np.where(times == times[0])[0][0]
        i_1 = np.where(times == times[-1])[0][0]
        intvl = [i_0, i_1]
   
    # V. sat completely out of times. Could be on either side
    else:
        intvl = None
    
    # intvl = interval
    return intvl


def read_aligned(ali_file=None):
    '''Read aligned data npz file'''
    paired_data = np.load(ali_file, allow_pickle=True)
    
    power   = paired_data['ref_p_aligned']
    times   = paired_data['time_array']

    return [power, times]


def noise_floor(sat_thresh, noi_thresh, data=None):
    '''Computes '''
    
    # compute the standard deviation of data, and use it to identify occupied channels
    σ = np.std(data)
    
    # Any channel with a max power >= σ has a satellite
    sat_cut = sat_thresh*σ
    chans_pow_max = np.amax(data, axis=0)
    
    # Exclude the channels with sats, to only have noise data
    noise_chans = np.where(chans_pow_max < sat_cut)[0]
    noise_data = data[:, noise_chans]
    
    # noise median, noise mad, noise threshold = μ + 3*σ
    μ_noise = np.median(noise_data)
    σ_noise = mad(noise_data, axis=None)
    # because we rescale power to have zero median
    noise_threshold = (μ_noise-μ_noise) + noi_thresh*σ_noise
    
    # scale the data so that it has zero median
    data = data - μ_noise
    
    return (data, noise_threshold)


def find_sat_channel(norad_id):
    '''All the magic [filtering] happens here'''

    # list of all occupied channels
    chans = []

    # for a sat window with multiple channels identified, use ephemeris to deterime which sats they could be

    # list of all sats identified
    sats = []

    # Lets make a counter of sats. If one or more channels are identified in a window, increment sat count by one

    sat_count = 0
    
   
    # Loop through days
    for day in range(len(dates)):
        
        # Loop through each 30 min obs in day
        for window in range(len(date_time[day])):
            
            ref_file = f'{ali_dir}/{dates[day]}/{date_time[day][window]}/rf0XX_S06XX_{date_time[day][window]}_aligned.npz'
            chrono_file = f'{chrono_dir}/{date_time[day][window]}.json'
            
            try:
                if Path(ref_file).is_file() and Path(chrono_file).is_file():
    
                    power, times = read_aligned(ali_file=ref_file)

                    # Scale noise floor to zero and determine noise threshold
                    power, noise_threshold = noise_floor(sat_thresh, noi_thresh, power)
                    
            except Exception:
                print(f'{date_time[day][window]}: ref file not found')
                continue

            try:    
                with open(chrono_file) as chrono:
                    chrono_ephem = json.load(chrono)
            
            except Exception:
                print(f'{date_time[day][window]}: chrono file not found')
                continue

            num_passes = len(chrono_ephem)
    
            norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(num_passes)]

            if norad_id in norad_list:
                norad_index = norad_list.index(norad_id)
            
                norad_ephem = chrono_ephem[norad_index]
                    
                rise_ephem  = norad_ephem["time_array"][0] 
                set_ephem   = norad_ephem["time_array"][-1]
                
                ## Max altitude that the sat attains
                #sat_alt_max = np.amax(norad_ephem["sat_alt"])
            
                ## Altitude threshold. Sats below 20 degrees are bound to be faint, and not register
                #if sat_alt_max >= alt_thresh:
                    
                Path(f'{plt_dir}/{norad_id}').mkdir(parents=True, exist_ok=True)
                   
                intvl = time_filter(rise_ephem, set_ephem, np.asarray(times))

                if intvl != None:
                    
                    w_start, w_stop = intvl
            
                    # length of sat pass. Only consider passes longer than 2 minutes
                    window_len = w_stop - w_start + 1
                    
                    # Slice [crop] the power/times arrays to the times of sat pass
                    power_c = power[w_start:w_stop+1, :]
                    times_c = times[w_start:w_stop+1]
                    
                    possible_chans = []

                    # Loop over every channel
                    for s_chan in range(len(power_c[0])):
                        
                        channel_power = power_c[:, s_chan]
                        
                        # Percentage of signal occupancy above noise threshold
                        window_occupancy = (np.where(channel_power >= noise_threshold))[0].size/window_len
            
                        max_s = np.amax(channel_power)
                        min_s = np.amin(channel_power)

                        # Arbitrary threshold below which satellites aren't counted
                        # Only continue if there is signal for more than 80% of satellite pass
                        if (max(channel_power) >= arb_thresh and 
                                occ_thresh <= window_occupancy < 1.00):
                         
                            ## fist and last 10 steps must be below the noise threshold
                            #if (all(p < noise_threshold for p in channel_power[:10]) and 
                            #    all(p < noise_threshold for p in channel_power[-11:-1])) is True:

                            #    # Center of gravity section
                            #    # Checks how central the signal is within the window
                            #    center, cog, frac_cen_offset = center_of_gravity(channel_power, times_c)

                            #    # Another threshold
                            #    # The Center of Gravity of signal is within 5% of center
                            #    if frac_cen_offset <= cog_thresh:
                            #    
                            #        # Plots the channel with satellite pass
                            #        plt_channel(
                            #                f'{plt_dir}/{norad_id}', times_c, power_c[:, s_chan],
                            #                s_chan, min_s, max_s, noise_threshold,
                            #                arb_thresh,center, cog, cog_thresh,
                            #                norad_id, f'{date_time[day][window]}')
                            #        
                            #        possible_chans.append(s_chan)
                            possible_chans.append(s_chan)
                    
                    # If channels are identified in the 30 min obs
                    n_chans = len(possible_chans)
                    
                    if n_chans > 0:

                        sat_count += 1
                        
                        # plot waterfall with sat window and all selected channels highlighted
                        plt_waterfall_pass(
                                f'{plt_dir}/{norad_id}', power, norad_id,
                                w_start, w_stop, possible_chans,
                                f'{date_time[day][window]}', cmap)
                        
                        # Add possible chans to ultimate list of chans, for histogram
                        chans.extend(possible_chans)

                        ## Plot ephemeris of lighly sats present in ephem window of norad_id
                        #plt_ids = []
                        #plt_alt = []
                        #plt_az  = []

                        #other_passes = []
                       
                        ## loop through all sats in chrono_ephem
                        #for s in range(len(chrono_ephem)):
                        #    times_sat = chrono_ephem[s]["time_array"]
                        #    
                        #    # Crop ephem of all sats to size of norad_id sat
                        #    intvl_ephem = time_filter(times_c[0], times_c[-1], np.asarray(times_sat))
                        #    
                        #    if intvl_ephem != None:
                        #        e_0, e_1 = intvl_ephem
                        #        
                        #        # altittude and occupancy filters
                        #        sat_alt_max = np.amax(chrono_ephem[s]["sat_alt"][e_0:e_1+1])
                        #        if sat_alt_max >= alt_thresh:
                        #           
                        #            if len(chrono_ephem[s]["sat_alt"][e_0:e_1+1]) >= occ_thresh*len(times_c):
                        #                
                        #                if chrono_ephem[s]["sat_id"][0] == norad_id:
                        #                
                        #                    plt_ids.extend(chrono_ephem[s]["sat_id"])
                        #                    plt_alt.append(chrono_ephem[s]["sat_alt"][e_0:e_1+1])
                        #                    plt_az.append(chrono_ephem[s]["sat_az"][e_0:e_1+1])
                        #                
                        #                else:
                        #                    other_ephem = []
                        #                    other_ephem.append(len(chrono_ephem[s]['sat_alt'][e_0:e_1+1]))
                        #                    other_ephem.extend(chrono_ephem[s]["sat_id"])
                        #                    other_ephem.append(chrono_ephem[s]["sat_alt"][e_0:e_1+1])
                        #                    other_ephem.append(chrono_ephem[s]["sat_az"][e_0:e_1+1])

                        #                    other_passes.append(other_ephem)


                        #other_passes = sorted(other_passes, key=lambda x: x[0])
                        #if n_chans > 1:
                        #    for e in other_passes[-(n_chans-1):][::-1]:   # BEWARE!!!! If two elements have same lenght, are they switched by reversing??
                        #        plt_ids.append(e[1])
                        #        plt_alt.append(e[2])
                        #        plt_az.append(e[3])
                        #
                        ## Plot sat ephemeris 
                        #sat_plot(f'{plt_dir}/{norad_id}', plt_ids, norad_id, plt_alt, plt_az, len(plt_ids), f'{date_time[day][window]}', 'passes')
                        sats.extend(plt_ids) 

    if chans != []:

        values, counts = np.unique(chans, return_counts=True)
        
        # if more than one channel has the same max number of passes, return all
        pop_chans = [int(values[i]) for i, j in enumerate(counts) if j == max(counts)]
        
        sat_data = {
                'sat_id': norad_id,
                'pop_chans': pop_chans,
                'sat_count': sat_count,
                'chans': chans,
                'sats': sats
                }
                
        Path(f'{out_dir}/channel_data/').mkdir(parents=True, exist_ok=True)

        with open(f'{out_dir}/channel_data/{int(norad_id)}.json','w') as f: 
            json.dump(sat_data, f, indent=4) 

if __name__=="__main__":

    import sys
    import argparse
    from pathlib import Path
    import concurrent.futures
    from channels_plt import plt_waterfall_pass, plt_channel, plt_hist, sat_plot
    
    
    parser = argparse.ArgumentParser(description="""
        Determines which channel each satellite occupies
        """)
    
    parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
    parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_channels/',help='Output directory. Default=./../../outputs/sat_channels/')
    parser.add_argument('--ali_dir', metavar='\b', default='./../../outputs/align_data/',help='Output directory. Default=./../../outputs/align_data/')
    parser.add_argument('--plt_dir', metavar='\b', default='./../../outputs/sat_channels/sat_passes',help='Output directory. Default=./../../outputs/sat_channels/sat_passes')
    parser.add_argument('--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
    parser.add_argument('--noi_thresh', metavar='\b', default=3,help='Noise Threshold: Multiples of MAD. Default=3.')
    parser.add_argument('--sat_thresh', metavar='\b', default=1,help='1 σ threshold to detect sats Default=1.')
    parser.add_argument('--arb_thresh', metavar='\b', default=12,help='Arbitrary Threshold to detect sats. Default=12 dB.')
    parser.add_argument('--alt_thresh', metavar='\b', default=20,help='Altitude Threshold to detect sats. Default=20 degrees.')
    parser.add_argument('--cog_thresh', metavar='\b', default=0.05,help='Center of Gravity Threshold to detect sats. Default=0.05')
    parser.add_argument('--occ_thresh', metavar='\b', default=0.80,help='Occupation Threshold of sat in window. Default=0.80')

    
    args = parser.parse_args()
    
    chrono_dir =        args.chrono_dir
    start_date =        args.start_date
    stop_date =         args.stop_date
    out_dir =           args.out_dir
    ali_dir =           args.ali_dir
    plt_dir =           args.plt_dir
    noi_thresh =        args.noi_thresh 
    sat_thresh =        args.sat_thresh 
    arb_thresh =        args.arb_thresh
    alt_thresh =        args.alt_thresh
    cog_thresh =        args.cog_thresh
    occ_thresh =        args.occ_thresh
    
    # Save logs 
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f'{out_dir}/logs_{start_date}_{stop_date}.txt', 'a')
    
    # Import list of tile names from rf_data.py
    tiles = rf.tile_names()
    
    
    # Import list of Norad catalogue IDs
    sat_list = [id for id in sat_ids.norad_ids.values()]
    
    # dates: list of days
    # date_time = list of 30 min observation windows
    dates, date_time = time_tree(start_date, stop_date)
        
    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(find_sat_channel, sat_list)
    
    for result in results:
        print(result)


