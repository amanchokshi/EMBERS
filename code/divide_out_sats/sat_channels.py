import math
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad

import sys
sys.path.append('../decode_rf_data')
from rf_data import read_data, plot_waterfall
from colormap import spectral


from scipy import interpolate
from scipy.signal import savgol_filter


parser = argparse.ArgumentParser(description="""
    Time alignes and smoothes all data 
    between two date ranges. Saves data
    pairs in organised directory 
    structure as .npz files.
    """)

#parser.add_argument('--data_dir', metavar='\b', help='Dir where date is saved')
#parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
#parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
#parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/align_data/',help='Output directory. Default=./../../outputs/align_data/')
parser.add_argument('--savgol_window', metavar='\b', default=151,help='Length of savgol window. Must be odd. Default=151')
parser.add_argument('--polyorder', metavar='\b', default=1,help='Order of polynomial to fit to savgol window. Default=1')
parser.add_argument('--interp_type', metavar='\b', default='cubic',help='Type of interpolation. Ex: cubic, linear, etc. Default=cubic')
parser.add_argument('--interp_freq', metavar='\b', default=2,help='Frequency at which to resample smoothed data, in Hertz. Default=2')


args = parser.parse_args()

#data_dir =          args.data_dir
#start_date =        args.start_date
#stop_date =         args.stop_date
#out_dir =           args.out_dir
savgol_window =     args.savgol_window
polyorder =         args.polyorder
interp_type =       args.interp_type
interp_freq =       args.interp_freq


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
    
    power, times = read_data(ref)
    savgol_pow = savgol_filter(power, savgol_window, polyorder, axis=0)
    time_smooth = np.arange(math.ceil(times[0]), math.floor(times[-1]), (1 / interp_freq))
    ref_sav_interp = interpolate.interp1d(times, savgol_pow, axis=0, kind=interp_type)
    power_smooth = ref_sav_interp(time_smooth)
    
    return (power_smooth, time_smooth)


def plt_waterfall_pass(power, sat_id, start, stop, idx):
    
    # Custom spectral colormap
    cmap = spectral()
    
    plt.style.use('dark_background')
    fig = plt.figure(figsize = (7,10))
    ax = fig.add_axes([0.12, 0.1, 0.72, 0.85])
    im = ax.imshow(power, interpolation='none', cmap=cmap)
    cax = fig.add_axes([0.88, 0.1, 0.03, 0.85])
    fig.colorbar(im, cax=cax)
    ax.set_aspect('auto')
    ax.set_title(f'Waterfall Plot: {sat_id}')
    
    ax.set_xlabel('Freq Channel')
    ax.set_ylabel('Time Step')
    
    ax.fill_between(range(len(power[0, :])),
            np.full(len(power[0, :]), start),
            np.full(len(power[0, :]), stop),
            color='white',
            alpha=0.2
            )
    
    #plt.show()
    plt.savefig(f'test/{idx}_{sat_id}.png')
    plt.close()


def plt_channel(times, channel_power, chan_num, sat_id, idx):
    plt.style.use('dark_background')
    plt.rcParams.update({"axes.facecolor": "#242a3c"})

    plt.scatter(times, channel_power, marker='.', alpha=0.2, color='#db3751', label='Data')
    plt.scatter(times[::49], np.full(len(times[::49]), noise_threshold),
            alpha=0.7, marker='.', color='#5cb7a9',
            label=f'Noise Cut: {noise_threshold:.2f} dBm')
    plt.scatter(times[::49], np.full(len(times[::49]), arbitrary_threshold),
            alpha=0.7, marker='.', color='#fba95f',
            label=f'Arbitrary Cut: {arbitrary_threshold} dBm')

    plt.ylim([min_s - 1, max_s + 1])
    plt.ylabel('Power [dBm]')
    plt.xlabel('Time [s]')
    plt.title(f'Satellite Pass in Channel: [{chan_num}]')
    plt.tight_layout()
    #leg = plt.legend(loc="upper right", frameon=True)
    leg = plt.legend(frameon=True)
    leg.get_frame().set_facecolor('white')
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)
    plt.savefig(f'test/{j}_{sat_id}_{chan_num}.png')
    plt.close()
    #plt.show()
    #break


ref_file = '../../../tiles_data/rf0XX/2019-10-07/rf0XX_2019-10-07-00:00.txt'
chrono_file = '../../outputs/sat_ephemeris/chrono_json/2019-10-07-00:00.json'


power, times = savgol_interp(ref_file, savgol_window, polyorder, interp_type, interp_freq )

# To first order, let us consider the median to be the noise floor
noise_f = np.median(power)
noise_mad = mad(power, axis=None)
noise_threshold = 3*noise_mad
arbitrary_threshold = 6 #dBm


# Scale the power to bring the median noise floor down to zero
power = power - noise_f
max_s = np.amax(power)
min_s = np.amin(power)


#plot_waterfall(power, times, 'rf0XX')
#plt.savefig(f'test/waterfall.png')
#plt.close()

used_channels = []

with open(chrono_file) as chrono:
    chrono_ephem = json.load(chrono)

    # Not sure how lambda functions work, but this does seem to sort chrono_ephem according to pass lenght
    chrono_ephem = sorted(chrono_ephem, key=lambda k: (k["time_array"][0] - k["time_array"][-1]))

    for j in range(len(chrono_ephem)):
    
        # Need to make sure that I'm not out of range on either side

        rise_ephem = chrono_ephem[j]["time_array"][0] 
        set_ephem  = chrono_ephem[j]["time_array"][-1]
        sat_id = chrono_ephem[j]["sat_id"][0]

        if sat_id == 44018:

            # Case I: Sat rises before obs starts and sets after the obs starts
            if rise_ephem <= times[0] and set_ephem >= times[0]:
                #print(f'{sat_id}: I')
                w_start = 0
                w_stop = list(times).index(set_ephem)

            # Case II: Sat rises after obs starts and sets before obs ends
            elif set_ephem < times[-1] and rise_ephem >= times[0]:
                #print(f'{sat_id}: II')
                w_start = list(times).index(rise_ephem)
                w_stop = list(times).index(set_ephem)

            # Case III: Sat rises before obs ends and sets after
            elif set_ephem >= times[-1] and rise_ephem <= times[-1]:
                #print(f'{sat_id}: III')
                w_start = list(times).index(rise_ephem)
                w_stop = list(times).index(times[-1])
            
            # Sat out of bounds
            else:
                #print(f'{sat_id}: out of bounds')
                w_start = 0
                w_stop = 0

            #plt_waterfall_pass(power, sat_id, w_start, w_stop, j)


            # length of sat pass
            window_len = w_stop - w_start + 1
            
            # Slice [crop] the power/times arrays to the times of sat pass
            power_c = power[w_start:w_stop+1, :]
            times_c = times[w_start:w_stop+1]

            occu_list = []
            chan_list = []
            
            for i in range(len(power_c[0])):
                
                if i not in used_channels:

                    channel_power = power_c[:, i]

                    # Arbitrary threshold below which satellites aren't counted
                    if max(channel_power) >= arbitrary_threshold:
                      
                        # Percentage of signal occupancy above noise threshold
                        window_occupancy = len([p for p in channel_power if p >= noise_threshold])/window_len
                        
                        # Only continue if there is signal for more than 80% of satellite pass
                        if window_occupancy >= 0.90:

                            occu_list.append(window_occupancy)
                            chan_list.append(i)


            if occu_list == []:
                pass
            else:
                # Channel number with max occupancy
                s_chan = chan_list[occu_list.index(max(occu_list))]
                
                # Exclude channels on either side of pass
                used_channels.extend([s_chan-1, s_chan, s_chan+1])
                print(f'Satellite {sat_id} in channel: {s_chan}, occupancy: {max(occu_list)*100:.2f}%')
                # Plots the channel with satellite pass
                #plt_channel(times_c, power_c[:, s_chan], s_chan, sat_id, j)

#TODO If the potential sat occupies more than 80% of the sat pass length, classify it as a sat!
#TODO Or, come up with alternative thresholding scheme.
#TODO Exclude that channel from next loop
#TODO Find a way to include more channels - weather sats

    
