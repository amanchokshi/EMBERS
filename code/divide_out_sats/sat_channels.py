import math
import json
import argparse
import numpy as np
import seaborn as sns
import concurrent.futures
#import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad
from datetime import datetime, timedelta
from pathlib import Path


import matplotlib
# Force matplotlib to not use X-Server backend
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
sys.path.append('../decode_rf_data')
sys.path.append('../sat_ephemeris')
import rf_data as rf
from colormap import spectral
import sat_ids


from scipy import interpolate
from scipy.signal import savgol_filter


parser = argparse.ArgumentParser(description="""
    Determines which channel each satellite occupies
    """)

parser.add_argument('--data_dir', metavar='\b', help='Dir where date is saved')
parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/divide_out_sats/sat_channels/',help='Output directory. Default=./../../outputs/divide_out_sats/sat_channels/')
parser.add_argument('--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
parser.add_argument('--savgol_window', metavar='\b', default=151,help='Length of savgol window. Must be odd. Default=151')
parser.add_argument('--polyorder', metavar='\b', default=1,help='Order of polynomial to fit to savgol window. Default=1')
parser.add_argument('--interp_type', metavar='\b', default='cubic',help='Type of interpolation. Ex: cubic, linear, etc. Default=cubic')
parser.add_argument('--interp_freq', metavar='\b', default=2,help='Frequency at which to resample smoothed data, in Hertz. Default=2')
parser.add_argument('--parallel', metavar='\b', default=True,help='If parallel=False, paralellization will be disabled.')



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


def plt_waterfall_pass(power, sat_id, start, stop, chs, date):
    '''Plot waterfall with sat window and occupied channels
    
    Args:
        power:          RF power array
        sat_id:         Norad cat ID
        start:          Start of epehm for sat_id
        stop:           Stop of ephem of sat_id
        chs:            Occupied channels [list]
        date:           Date of observation
    '''

    # Custom spectral colormap
    cmap = spectral()
    
    plt.style.use('dark_background')
    fig = plt.figure(figsize = (7,10))
    ax = fig.add_axes([0.12, 0.1, 0.72, 0.85])
    im = ax.imshow(power, interpolation='none', cmap=cmap)
    cax = fig.add_axes([0.88, 0.1, 0.03, 0.85])
    fig.colorbar(im, cax=cax)
    ax.set_aspect('auto')
    ax.set_title(f'Waterfall Plot: {sat_id} in {chs}')
    
    ax.set_xlabel('Freq Channel')
    ax.set_ylabel('Time Step')
    
    # horizontal highlight of ephem
    ax.axhspan(start, stop, alpha=0.1, color='white')
    
    # vertical highlight of channel
    for ch in chs:
        ax.axvspan(ch-1.0, ch+0.6, alpha=0.2, color='white')
    
    #plt.show()
    plt.savefig(f'{out_dir}/{sat_id}/{date}_{sat_id}_waterfall.png')
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)


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


def plt_channel(times, channel_power, chan_num, min_s, max_s, noise_threshold, arbitrary_threshold, center, cog, sat_id, date):
    '''Plot power in channel, with various thresholds
    
    Args:
        times:          Time array
        channel_power:  Power in channel
        chan_num:       Channel Number
        min_s:          Minimum signal in channel_power
        max_s:          Maximum signal in channel_power
        noise_threshold: Noise Threshold (n*MAD)
        arbitrary_threshold: Arbitrary threshold used to only select bright passes
        center:         Center of channel_power
        cog:            Center of gravity of channel_power
        sat_id:         Norad Cat ID
        date:           Date of observation
        '''
    
    
    plt.style.use('seaborn')
    
    # plt channel power
    plt.plot(times, channel_power, linestyle='-', linewidth=2, alpha=1.0, color='#db3751', label='Data')
    plt.fill_between(times, channel_power, color='#db3751', alpha=0.7)
    
    plt.axhline(arbitrary_threshold,alpha=1.0, linestyle='-', linewidth=2,
            color='#fba95f', label=f'Arbitrary Cut: {arbitrary_threshold} dBm')
    plt.axhspan(-1, arbitrary_threshold, color='#fba95f', alpha=0.4)
    
    plt.axhline(noise_threshold, linestyle='-', linewidth=2, color='#5cb7a9',
            label=f'Noise Cut: {noise_threshold:.2f} dBm')
    plt.axhspan(-1, noise_threshold, color='#5cb7a9', alpha=0.4)
    
    plt.axvspan(center-(0.02*len(times)), center+(0.02*len(times)), color='#2b2e4a', alpha=0.4)
    plt.axvline(center, color='#2b2e4a', alpha=1, label='Center ± 2% ')
    plt.axvline(cog, color='#8cba51', alpha=1, label='CoG')
    
    
    plt.ylim([min_s - 1, max_s + 1])
    plt.xlim([times[0], times[-1]])
    plt.ylabel('Power [dBm]')
    plt.xlabel('Time [s]')
    plt.title(f'Satellite Pass in Channel: [{chan_num}]')
    plt.tight_layout()
    leg = plt.legend(frameon=True)
    leg.get_frame().set_facecolor('grey')
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)
    plt.savefig(f'{out_dir}/{sat_id}/{date}_{sat_id}_{chan_num}_channel.png')
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)


def sat_plot(ids, norad_id, alt, az, num_passes, date):
    '''Plots satellite passes
    
    Args:
        alt: list of altitude values
        az: list of azimuth values
        num_passes: Number of satellite passes
    '''
    
    
    # Set up the polar plot.
    plt.style.use('dark_background')
    figure = plt.figure(figsize=(8,6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,30,60,90], angle=22)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title(f'Satellite on {date}: {num_passes} Passes', y=1.08)
    ax.grid(color='#8bbabb', linewidth=1.6, alpha=0.6)
    plt.tight_layout()
    
    for i in range(len(alt)):
        plt.plot(az[i], alt[i], '-', linewidth=1.6, label=f'{ids[i]}')
        plt.legend(bbox_to_anchor=(0.28, 1.0, 1., .102), loc='upper right')
    
    plt.savefig(f'{out_dir}/{norad_id}/{date}_passes.png')
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)


def plt_hist(chans, norad_id, pop_chan):
    '''Plt histogram of occupied channels'''
    
    #plt.rcParams.update(plt.rcParamsDefault)
    sns.set()
    values, counts = np.unique(chans, return_counts=True)
    y_pos = np.arange(len(values))
    plt.bar(y_pos, counts, color=sns.color_palette("GnBu_d", len(counts)))
    plt.ylabel('Number of Passes in Channel')
    plt.xlabel('Channel')
    plt.xticks(y_pos, values)
    plt.title(f'Possible Transmission Channels of Satellite [{norad_id}]')
    plt.tight_layout()
    plt.savefig(f'{out_dir}/{norad_id}/Channels_histogram_{norad_id}_{pop_chan}.png')
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)



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


                
                    for j in range(num_passes):
                    
                        rise_ephem  = chrono_ephem[j]["time_array"][0] 
                        set_ephem   = chrono_ephem[j]["time_array"][-1]
                        sat_alt_max = np.amax(chrono_ephem[j]["sat_alt"])
                        sat_id      = chrono_ephem[j]["sat_id"][0]
                
                        # Focus on one sat at a time
                        if sat_id == norad_id:
                            
                            # Altitude threshold. Sats below 20 degrees are bound to be faint, and not register
                            if sat_alt_max >= 20:
                                
                                Path(f'{out_dir}/{sat_id}').mkdir(parents=True, exist_ok=True)
                                    
                                # window start, stop
                                # Case I: Sat rises after obs starts and sets before obs ends
                                if set_ephem <= times[-1] and rise_ephem >= times[0]:
                                    w_start = list(times).index(rise_ephem)
                                    w_stop = list(times).index(set_ephem)
                
                                # Case II: Sat rises before obs starts and sets after the obs starts
                                elif rise_ephem < times[0] and set_ephem > times[0]:
                                    # w_start = 0
                                    w_start = list(times).index(times[0])
                                    w_stop = list(times).index(set_ephem)
                
                                # Case III: Sat rises before obs ends and sets after
                                elif set_ephem > times[-1] and rise_ephem < times[-1]:
                                    w_start = list(times).index(rise_ephem)
                                    w_stop = list(times).index(times[-1])
                                
                                # Sat out of bounds
                                else:
                                    w_start = 0
                                    w_stop = 0
                
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

                                                    # Another threshold
                                                    # The Center of Gravity of signal is within 2% of center
                                                    if frac_cen_offset <= 0.02:
                                                    
                                                        # Plots the channel with satellite pass
                                                        plt_channel(times_c, power_c[:, s_chan],
                                                                s_chan, min_s, max_s, noise_threshold,
                                                                arbitrary_threshold,center, cog, sat_id, f'{date_time[day][window]}')
                                                        
                                                        possible_chans.append(s_chan)
                                                        #chans.append(s_chan)
                                    
                                    if len(possible_chans) < 1:
                                        pass
                                    else:
                                        plt_waterfall_pass(power, sat_id, w_start, w_stop, possible_chans, f'{date_time[day][window]}')
                                        chans.extend(possible_chans)

                                        wt_start = times_c[0]
                                        wt_stop = times_c[-1]

                                        # TODO only want to plot ephem of sats within the window of w_start:w_stop+1
                                        
                                        ids = [chrono_ephem[i]["sat_id"][0] for i in range(num_passes)]
                                        alt = [chrono_ephem[i]["sat_alt"] for i in range(num_passes)]
                                        az  = [chrono_ephem[i]["sat_az"] for i in range(num_passes)]

                                        sat_plot(ids, norad_id, alt, az, num_passes, f'{date_time[day][window]}')

                                    

    if chans == []:
        print(f'No identified sat channels for sat  {norad_id}')
    else:
        chans = np.array(chans).astype('int64')
        counts = np.bincount(chans)
        pop_chan = np.argmax(counts)
        plt_hist(chans, norad_id, pop_chan)
        print(f'Most frequently occupied channel for sat {norad_id}: {pop_chan}')


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


# Time stuff to help traverse data dir tree.
t_start = datetime.strptime(start_date, '%Y-%m-%d')
t_stop = datetime.strptime(stop_date, '%Y-%m-%d')
n_days = (t_stop - t_start).days

dates = []
date_time = []

for i in range(n_days+1):
    day = t_start + timedelta(days=i)
    date = day.strftime('%Y-%m-%d')
    dates.append(date)
    d_t = []
    
    for j in range(48):
        t_delta = datetime.strptime(date,'%Y-%m-%d') + timedelta(minutes=30*j)
        d_time = t_delta.strftime('%Y-%m-%d-%H:%M')
        d_t.append(d_time)

    date_time.append(d_t)    


data_dir = Path(data_dir)
chrono_dir = Path(chrono_dir)
out_dir = Path(out_dir)
    
#norad_id = 41180
#find_sat_channel(norad_id)

if parallel != True:
    for norad_id in sat_list:
        find_sat_channel(norad_id)
        break
else:
    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor(max_workers=40) as executor:
        results = executor.map(find_sat_channel, sat_list)

    for result in results:
        print(result)


