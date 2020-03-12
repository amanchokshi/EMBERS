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
            
            ref_file = f'{data_dir}/{ref_tile}/{dates[day]}/{ref_tile}_{date_time[day][window]}.txt'
            chrono_file = f'{chrono_dir}/{date_time[day][window]}.json'
            
            try:
                if Path(ref_file).is_file() and Path(chrono_file).is_file():
    
                    power, times = savgol_interp(ref_file, savgol_window, polyorder, interp_type, interp_freq )
                    
                    # To first order, let us consider the median to be the noise floor
                    noise_f = np.median(power)
                    noise_mad = mad(power, axis=None)
                    noise_threshold = noi_thresh*noise_mad
                    arbitrary_threshold = 10 #dBm
                    
                    # Scale the power to bring the median noise floor down to zero
                    power = power - noise_f
            
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
                
                # Max altitude that the sat attains
                sat_alt_max = np.amax(norad_ephem["sat_alt"])
            
                # Altitude threshold. Sats below 20 degrees are bound to be faint, and not register
                if sat_alt_max >= alt_thresh:
                    
                    Path(f'{out_dir}/{norad_id}').mkdir(parents=True, exist_ok=True)
                       
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
                             
                                # fist and last 10 steps must be below the noise threshold
                                if (all(p < noise_threshold for p in channel_power[:10]) and 
                                    all(p < noise_threshold for p in channel_power[-11:-1])) is True:

                                    # Center of gravity section
                                    # Checks how central the signal is within the window
                                    center, cog, frac_cen_offset = center_of_gravity(channel_power, times_c)

                                    # Another threshold
                                    # The Center of Gravity of signal is within 5% of center
                                    if frac_cen_offset <= cog_thresh:
                                    
                                        # Plots the channel with satellite pass
                                        plt_channel(
                                                out_dir, times_c, power_c[:, s_chan],
                                                s_chan, min_s, max_s, noise_threshold,
                                                arbitrary_threshold,center, cog, cog_thresh,
                                                norad_id, f'{date_time[day][window]}')
                                        
                                        possible_chans.append(s_chan)
                        
                        # If channels are identified in the 30 min obs
                        n_chans = len(possible_chans)
                        
                        if n_chans > 0:

                            sat_count += 1
                            
                            # plot waterfall with sat window and all selected channels highlighted
                            plt_waterfall_pass(
                                    out_dir, power, norad_id,
                                    w_start, w_stop, possible_chans,
                                    f'{date_time[day][window]}', cmap)
                            
                            # Add possible chans to ultimate list of chans, for histogram
                            chans.extend(possible_chans)

                            # Plot ephemeris of lighly sats present in ephem window of norad_id
                            plt_ids = []
                            plt_alt = []
                            plt_az  = []

                            other_passes = []
                           
                            # loop through all sats in chrono_ephem
                            for s in range(len(chrono_ephem)):
                                times_sat = chrono_ephem[s]["time_array"]
                                
                                # Crop ephem of all sats to size of norad_id sat
                                intvl_ephem = time_filter(times_c[0], times_c[-1], np.asarray(times_sat))
                                
                                if intvl_ephem != None:
                                    e_0, e_1 = intvl_ephem
                                    
                                    # altittude and occupancy filters
                                    sat_alt_max = np.amax(chrono_ephem[s]["sat_alt"][e_0:e_1+1])
                                    if sat_alt_max >= alt_thresh:
                                       
                                        if len(chrono_ephem[s]["sat_alt"][e_0:e_1+1]) >= occ_thresh*len(times_c):
                                            
                                            if chrono_ephem[s]["sat_id"][0] == norad_id:
                                            
                                                plt_ids.extend(chrono_ephem[s]["sat_id"])
                                                plt_alt.append(chrono_ephem[s]["sat_alt"][e_0:e_1+1])
                                                plt_az.append(chrono_ephem[s]["sat_az"][e_0:e_1+1])
                                            
                                            else:
                                                other_ephem = []
                                                other_ephem.append(len(chrono_ephem[s]['sat_alt'][e_0:e_1+1]))
                                                other_ephem.extend(chrono_ephem[s]["sat_id"])
                                                other_ephem.append(chrono_ephem[s]["sat_alt"][e_0:e_1+1])
                                                other_ephem.append(chrono_ephem[s]["sat_az"][e_0:e_1+1])

                                                other_passes.append(other_ephem)


                            other_passes = sorted(other_passes, key=lambda x: x[0])
                            if n_chans > 1:
                                for e in other_passes[-(n_chans-1):][::-1]:   # BEWARE!!!! If two elements have same lenght, are they switched by reversing??
                                    plt_ids.append(e[1])
                                    plt_alt.append(e[2])
                                    plt_az.append(e[3])
                            
                            # Plot sat ephemeris 
                            sat_plot(out_dir, plt_ids, norad_id, plt_alt, plt_az, len(plt_ids), f'{date_time[day][window]}', 'passes')
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
    
    parser.add_argument('--data_dir', metavar='\b', help='Dir where date is saved')
    parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
    parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_channels/',help='Output directory. Default=./../../outputs/sat_channels/')
    parser.add_argument('--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
    parser.add_argument('--savgol_window', metavar='\b', default=151,help='Length of savgol window. Must be odd. Default=151')
    parser.add_argument('--polyorder', metavar='\b', default=1,help='Order of polynomial to fit to savgol window. Default=1')
    parser.add_argument('--interp_type', metavar='\b', default='cubic',help='Type of interpolation. Ex: cubic, linear, etc. Default=cubic')
    parser.add_argument('--interp_freq', metavar='\b', default=1,help='Frequency at which to resample smoothed data, in Hertz. Default=2')
    parser.add_argument('--parallel', metavar='\b', default=True,help='If parallel=False, paralellization will be disabled.')
    
    parser.add_argument('--noi_thresh', metavar='\b', default=10,help='Noise Threshold: Multiples of MAD. Default=10.')
    parser.add_argument('--arb_thresh', metavar='\b', default=12,help='Arbitrary Threshold to detect sats. Default=12 dB.')
    parser.add_argument('--alt_thresh', metavar='\b', default=20,help='Altitude Threshold to detect sats. Default=20 degrees.')
    parser.add_argument('--cog_thresh', metavar='\b', default=0.05,help='Center of Gravity Threshold to detect sats. Default=0.05')
    parser.add_argument('--occ_thresh', metavar='\b', default=0.80,help='Occupation Threshold of sat in window. Default=0.80')

    
    args = parser.parse_args()
    
    data_dir =          args.data_dir
    chrono_dir =        args.chrono_dir
    start_date =        args.start_date
    stop_date =         args.stop_date
    out_dir =           args.out_dir
    savgol_window =     args.savgol_window
    polyorder =         args.polyorder
    interp_type =       args.interp_type
    interp_freq =       args.interp_freq
    noi_thresh =        args.noi_thresh 
    arb_thresh =        args.arb_thresh
    alt_thresh =        args.alt_thresh
    cog_thresh =        args.cog_thresh
    occ_thresh =        args.occ_thresh
    parallel=           args.parallel
    
    # Save logs 
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f'{out_dir}/logs_{start_date}_{stop_date}.txt', 'a')
    
    # Import list of tile names from rf_data.py
    tiles = rf.tile_names()
    
    # Only using rf0XX for now
    ref_tile = tiles[0]
    
    # Import list of Norad catalogue IDs
    sat_list = [id for id in sat_ids.norad_ids.values()]
    
    
    # Path to important locations
    data_dir = Path(data_dir)
    chrono_dir = Path(chrono_dir)
    out_dir = Path(out_dir)
    
    # dates: list of days
    # date_time = list of 30 min observation windows
    dates, date_time = time_tree(start_date, stop_date)
        
    #norad_id = 41180
    #find_sat_channel(norad_id)
    
    if parallel != True:
        for norad_id in sat_list:
            find_sat_channel(norad_id)
            #break
    
    else:
        try:
        # Parallization magic happens here
            with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
                results = executor.map(find_sat_channel, sat_list)
        
            for result in results:
                print(result)
        except Exception:
            print('Inssuficient Memory. Try using --parallel=False flag to run code serially')


