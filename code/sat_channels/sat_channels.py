

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
   
    chans = []
   
    # Loop through days
    for day in range(len(dates)):
        
        # Loop through each 30 min obs in day
        for window in range(len(date_time[day])):
            
            ref_file = f'{data_dir}/{ref_tile}/{dates[day]}/{ref_tile}_{date_time[day][window]}.txt'
            chrono_file = f'{chrono_dir}/{date_time[day][window]}.json'
   
            if Path(ref_file).is_file() and Path(chrono_file).is_file():
    
                power, times = savgol_interp(ref_file, savgol_window, polyorder, interp_type, interp_freq )
                
                # To first order, let us consider the median to be the noise floor
                noise_f = np.median(power)
                noise_mad = mad(power, axis=None)
                noise_threshold = 10*noise_mad
                arbitrary_threshold = 10 #dBm
                
                
                # Scale the power to bring the median noise floor down to zero
                power = power - noise_f
                
                
                with open(chrono_file) as chrono:
                    chrono_ephem = json.load(chrono)

                    num_passes = len(chrono_ephem)
    
                    norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(num_passes)]

                    if norad_id in norad_list:
                        norad_index = norad_list.index(norad_id)
                    
                        norad_ephem = chrono_ephem[norad_index]
                            
                        rise_ephem  = norad_ephem["time_array"][0] 
                        set_ephem   = norad_ephem["time_array"][-1]
                        sat_id      = norad_ephem["sat_id"][0]
                        
                        # Max altitude that the sat attains
                        sat_alt_max = np.amax(norad_ephem["sat_alt"])
                
                        ## Focus on one sat at a time
                        #if sat_id == norad_id:
                            
                        # Altitude threshold. Sats below 20 degrees are bound to be faint, and not register
                        if sat_alt_max >= 20:
                            
                            Path(f'{out_dir}/{norad_id}').mkdir(parents=True, exist_ok=True)
                               
                            intvl = time_filter(rise_ephem, set_ephem, np.asarray(times))

                            if intvl != None:
                                
                                w_start, w_stop = intvl
                
                                # length of sat pass. Only consider passes longer than 2 minutes
                                window_len = w_stop - w_start + 1
                                
                                if window_len >= 240: #(240*0.5 = 120s)
    
                                    # Slice [crop] the power/times arrays to the times of sat pass
                                    power_c = power[w_start:w_stop+1, :]
                                    times_c = times[w_start:w_stop+1]
                                   
                                    possible_chans = []

                                    # Loop over every channel
                                    for s_chan in range(len(power_c[0])):
                                        
                                        channel_power = power_c[:, s_chan]
                
                                        max_s = np.amax(channel_power)
                                        min_s = np.amin(channel_power)

                                        # Arbitrary threshold below which satellites aren't counted
                                        if max(channel_power) >= arbitrary_threshold:
                                          
                                            # Percentage of signal occupancy above noise threshold
                                            window_occupancy = (np.where(channel_power >= noise_threshold))[0].size/window_len
                                            
                                            # Only continue if there is signal for more than 80% of satellite pass
                                            if window_occupancy >= 0.80 and window_occupancy < 1.00:

                                                # Make sure that the ends are close to the noise floor
                                                if (all(p < noise_threshold for p in channel_power[:10]) and 
                                                    all(p < noise_threshold for p in channel_power[-11:-1])) is True:

                                                    center, cog, frac_cen_offset = center_of_gravity(channel_power, times_c)

                                                    c_thresh = 0.05
                                                    # Another threshold
                                                    # The Center of Gravity of signal is within 5% of center
                                                    if frac_cen_offset <= c_thresh:
                                                    
                                                        # Plots the channel with satellite pass
                                                        plt_channel(
                                                                out_dir, times_c, power_c[:, s_chan],
                                                                s_chan, min_s, max_s, noise_threshold,
                                                                arbitrary_threshold,center, cog, c_thresh,
                                                                norad_id, f'{date_time[day][window]}')
                                                        
                                                        possible_chans.append(s_chan)
                                    
                                    # If channels are identified in the 30 min obs
                                    if len(possible_chans) > 0:

                                        plt_waterfall_pass(
                                                out_dir, power, norad_id,
                                                w_start, w_stop, possible_chans,
                                                f'{date_time[day][window]}', cmap)
                                        
                                        # Add possible chans to ultimate list of chans, for histogram
                                        chans.extend(possible_chans)


                                        # Plot ephemeris of all sats present in ephem window of norad_id sat
                                        ids = []
                                        alt = []
                                        az  = []
                                       
                                        # loop through all sats in chrono_ephem
                                        for s in range(len(chrono_ephem)):
                                            times_sat = chrono_ephem[s]["time_array"]
                                            
                                            # Crop ephem of all sats to size of norad_id sat
                                            intvl_ephem = time_filter(times_sat[0], times_sat[-1], times_c)
                                            
                                            if intvl_ephem != None:
                                                e_0, e_1 = intvl_ephem
                                                
                                                ids.extend(chrono_ephem[s]["sat_id"])
                                                alt.append(chrono_ephem[s]["sat_alt"][e_0:e_1+1])
                                                az.append(chrono_ephem[s]["sat_az"][e_0:e_1+1])
                                        
                                        if norad_id in ids:
                                            sat_plot(out_dir, ids, norad_id, alt, az, len(ids), f'{date_time[day][window]}')

                                    

    if chans == []:
        print(f'No identified sat channels for sat {norad_id}')
    else:
        chans = np.array(chans).astype('int64')
        counts = np.bincount(chans)
        pop_chan = np.argmax(counts)
        plt_hist(out_dir, chans, norad_id, pop_chan)
        print(f'Most frequently occupied channel for sat {norad_id}: {pop_chan}')


if __name__=="__main__":

    import math
    # Force matplotlib to not use X-Server backend
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import json
    import argparse
    import numpy as np
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
    import rf_data as rf
    from colormap import spectral
    # Custom spectral colormap
    cmap = spectral()
    import sat_ids
    
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
        
    norad_id = 41180
    find_sat_channel(norad_id)
    
    #if parallel != True:
    #    for norad_id in sat_list:
    #        find_sat_channel(norad_id)
    #        break
    #
    #else:
    #    # Parallization magic happens here
    #    with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
    #        results = executor.map(find_sat_channel, sat_list)
    #
    #    for result in results:
    #        print(result)


