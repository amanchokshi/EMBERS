import sys
import math
import json
import numpy as np
import seaborn as sns
from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.stats import median_absolute_deviation as mad
from datetime import datetime, timedelta
from pathlib import Path

# Force matplotlib to not use X-Server backend
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pl

## Custom spectral colormap
sys.path.append('../decode_rf_data')
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

def plt_waterfall_pass(out_dir, power, sat_id, start, stop, chs, date, cmap):
    '''Plot waterfall with sat window and occupied channels
    
    Args:
        power:          RF power array
        sat_id:         Norad cat ID
        start:          Start of epehm for sat_id
        stop:           Stop of ephem of sat_id
        chs:            Occupied channels [list]
        date:           Date of observation
    '''

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


def plt_channel(
        out_dir, times, channel_power,
        chan_num, min_s, max_s, noise_threshold,
        arbitrary_threshold, center, cog, c_thresh, sat_id, date):

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
        c_thresh:       Threshold about center
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
    
    plt.axvspan(center-(c_thresh*len(times)), center+(c_thresh*len(times)), color='#2b2e4a', alpha=0.4)
    plt.axvline(center, color='#2b2e4a', alpha=1, label=f'Center ± {c_thresh*100}% ')
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


def sat_plot(out_dir, ids, norad_id, alt, az, num_passes, date, name):
    '''Plots satellite passes
    
    Args:
        alt: list of altitude values
        az: list of azimuth values
        num_passes: Number of satellite passes
    '''
    
    
    # Set up the polar plot.
    plt.style.use('seaborn')
    plt.style.use('dark_background')
    figure = plt.figure(figsize=(7,6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,30,60,90], angle=22)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title(f'{num_passes} satellite passes in {date} [{norad_id}] window', y=1.05)
    ax.grid(color='grey', linewidth=1.6, alpha=0.5)

    colors = pl.cm.Spectral(np.linspace(0.17,0.9,len(alt)))
    
    for i in range(len(alt)):
        plt.plot(az[i], alt[i], '-', linewidth=2.4, alpha=0.8, color=colors[i], label=f'{ids[i]}')
        plt.legend()
    
    leg = plt.legend(frameon=True, bbox_to_anchor=(0.28, 1.0, 1., -0.95), loc='center right', title="Norad SatID")
    leg.get_frame().set_facecolor('grey')
    leg.get_frame().set_alpha(0.4)
    for l in leg.legendHandles:
        l.set_alpha(1)
    
    plt.tight_layout()
    plt.savefig(f'{out_dir}/{norad_id}/{date}_{norad_id}_{name}.png')
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)


