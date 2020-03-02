import math
import json
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad

import sys
sys.path.append('../decode_rf_data')
from rf_data import read_data, plot_waterfall


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


def plt_channel(times, channel_power, chan_num):
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
    plt.savefig(f'test/channel_{chan_num}.png')
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


plot_waterfall(power, times, 'rf0XX')
plt.savefig(f'test/waterfall.png')
plt.close()

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
        print(sat_id)

        if rise_ephem < times[0]:
            print('I')
            start_ind = 0
            print(times[0])
            stop_ind = list(times).index(set_ephem)

        elif set_ephem > times[-1]:
            print('II')
            start_ind = list(times).index(rise_ephem)
            stop_ind = -2

        else:
            print('III')
            start_ind = list(times).index(rise_ephem)
            stop_ind = list(times).index(set_ephem)
        
        # length of sat pass
        window_len = stop_ind - start_ind
        
        # Slice the power/times arrays to the times of sat pass
        power = power[start_ind:stop_ind+1, :]
        times = times[start_ind:stop_ind+1]

        occu_list = []
        chan_list = []
        
        for i in range(len(power[0])):
            
            if i not in used_channels:

                channel_power = power[:, i]

                # Arbitrary threshold below which satellites aren't counted
                if max(channel_power) >= arbitrary_threshold:
                  
                    # Percentage of signal occupancy above noise threshold
                    window_occupancy = len([p for p in channel_power if p >= noise_threshold])/window_len
                    
                    # Only continue if there is signal for more than 80% of satellite pass
                    if window_occupancy >= 0.8:

                        occu_list.append(window_occupancy)
                        chan_list.append(i)


        if occu_list == []:
            pass
        else:
            print(occu_list)
            # Channel number with max occupancy
            s_chan = chan_list[occu_list.index(max(occu_list))]
            
            # Exclude channels on either side of pass
            used_channels.extend([s_chan-1, s_chan, s_chan+1])
            print(f'Satellite in channel: {s_chan}, occupancy: {max(occu_list)*100:.2f}%')

            # Plots the channel with satellite pass
            plt_channel(times, power[:, s_chan], s_chan)


#TODO If the potential sat occupies more than 80% of the sat pass length, classify it as a sat!
#TODO Or, come up with alternative thresholding scheme.
#TODO Exclude that channel from next loop
#TODO Find a way to include more channels - weather sats

    
