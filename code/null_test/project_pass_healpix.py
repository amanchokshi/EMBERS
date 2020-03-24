import sys
import json
import argparse
import numpy as np
import healpy as hp
from pathlib import Path
import concurrent.futures
import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad

sys.path.append('../sat_channels')
from sat_channels import time_tree, time_filter

# Custom spectral colormap
sys.path.append('../decode_rf_data')
from colormap import spectral
cmap = spectral()

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


def plt_channel(
        out_dir, times, channel_power,
        chan_num, min_s, max_s, noise_threshold,
        sat_id, date):

    '''Plot power in channel, with various thresholds
    
    Args:
        times:          Time array
        channel_power:  Power in channel
        chan_num:       Channel Number
        min_s:          Minimum signal in channel_power
        max_s:          Maximum signal in channel_power
        noise_threshold: Noise Threshold (n*MAD)
        sat_id:         Norad Cat ID
        date:           Date of observation
        '''
    
    
    plt.style.use('seaborn')
    
    # plt channel power
    plt.plot(times, channel_power, linestyle='-', linewidth=2, alpha=1.0, color='#729d39', label='Data')
    plt.fill_between(times, channel_power, color='#83b271', alpha=0.7)
    
    plt.axhline(noise_threshold, linestyle='-', linewidth=2, color='#ff6138',
            label=f'Noise Cut: {noise_threshold:.2f} dBm')
    plt.axhspan(-1, noise_threshold, color='#ff895d', alpha=0.4)
    
    plt.ylim([min_s - 1, max_s + 1])
    plt.xlim([times[0], times[-1]])
    plt.ylabel('Power [dBm]')
    plt.xlabel('Time [s]')
    plt.title(f'Satellite [{sat_id}] Pass @ {date} in Channel: [{chan_num}]')
    plt.tight_layout()
    leg = plt.legend(frameon=True)
    leg.get_frame().set_facecolor('grey')
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)
    plt.savefig(f'{out_dir}/{date}_{sat_id}_{chan_num}_channel.png')
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)



def power_ephem(
        ref,
        ref_file,
        chrono_file,
        sat_id,
        chan,
        timestamp
        ):

    '''Create power, alt, az arrays at constant cadence'''

    power, times = read_aligned(ali_file=ref_file)

    # Scale noise floor to zero and determine noise threshold
    power, noise_threshold = noise_floor(sat_thresh, noi_thresh, power)

    with open(chrono_file) as chrono:
        chrono_ephem = json.load(chrono)
    
        norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(len(chrono_ephem))]
        
        norad_index = norad_list.index(sat_id)

        norad_ephem = chrono_ephem[norad_index]
            
        rise_ephem  = norad_ephem["time_array"][0] 
        set_ephem   = norad_ephem["time_array"][-1]
        
        intvl = time_filter(rise_ephem, set_ephem, np.asarray(times))
        
        if intvl != None:
            
            w_start, w_stop = intvl
            
            # Slice [crop] the power/times arrays to the times of sat pass
            channel_power = power[w_start:w_stop+1, chan]
            times_c = times[w_start:w_stop+1]

            alt = np.asarray(norad_ephem["sat_alt"])
            az  = np.asarray(norad_ephem["sat_az"])

            if np.where(channel_power >= noise_threshold)[0].size != 0: 
                good_power = channel_power[np.where(channel_power >= noise_threshold)[0]]
                good_alt = alt[np.where(channel_power >= noise_threshold)[0]]
                good_az  = az[np.where(channel_power >= noise_threshold)[0]]
                
                
                if plots == 'True':
                    plt_channel(f'{plt_dir}/{ref}', times_c, channel_power, chan, min(channel_power), max(channel_power), noise_threshold, sat_id, timestamp )
                
                return [good_power, good_alt, good_az]

            else:
                return 0

        else:
            return 0
                

def proj_ref_healpix(ref):
    
    if plots == 'True':
        Path(f'{plt_dir}/{ref}').mkdir(parents=True, exist_ok=True)

    # Initialize empty beam map and list
    # a list of empty lists to append values to
    ref_map  = [[] for pixel in range(hp.nside2npix(nside))]

    # a list of zeros, to increment when a new values is added
    ref_counter=np.zeros(hp.nside2npix(nside))
    
    # Help traverse all 30 min obs b/w start & stop
    dates, date_time = time_tree(start_date, stop_date)

    # Loop through days
    for day in range(len(dates)):
        
        # Loop through each 30 min obs in day
        for window in range(len(date_time[day])):

            date = dates[day]
            timestamp = date_time[day][window]
            
            if 'XX' in ref:
                ref_file = f'{ali_dir}/{date}/{timestamp}/{ref}_S06XX_{timestamp}_aligned.npz'
            else:
                ref_file = f'{ali_dir}/{date}/{timestamp}/{ref}_S06YY_{timestamp}_aligned.npz'
            
            
            chrono_file = f'{chrono_dir}/{timestamp}.json'
            channel_map = f'{map_dir}/{timestamp}.json'

            
            try:
                Path(ref_file).is_file()

                with open(chrono_file) as chrono:
                    chrono_ephem = json.load(chrono)

                    if chrono_ephem != []:
                
                        norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(len(chrono_ephem))]

                        if norad_list != []:

                            with open(channel_map) as ch_map:
                                chan_map = json.load(ch_map)

                                chan_sat_ids = [int(i) for i in list(chan_map.keys())]
                        
                                for sat in chan_sat_ids:

                                    if sat in norad_list:

                                        chan = chan_map[f'{sat}']

                                        sat_data = power_ephem(
                                                ref,
                                                ref_file,
                                                chrono_file,
                                                sat,
                                                chan,
                                                timestamp
                                                )
                                        
                                        if sat_data != 0:
                                            channel_power, alt, az = sat_data

                                            # Altitude is in deg while az is in radians
                                            # convert alt to radians
                                            alt = np.radians(alt)
                                            az  = np.asarray(az)

                                            # To convert from Alt/Az to θ/ϕ spherical coordinates
                                            # Jack's convention, not sure about ɸ
                                            # θ = 90 - Alt
                                            # ɸ = 180 - Az

                                            # Healpix uses sperical coordinates
                                            θ = np.pi/2 - alt
                                            ɸ = np.pi - az

                                            # Since we need to slice along NS & EW, and nside = 32 healpix does not 
                                            # straight lines of pixels vertically or horizontally, but it does have
                                            # them diagonally. We rotate ɸ by 45° to be able to slice NS & EW
                                            ɸ_rot = ɸ + (np.pi / 4)

                                            # Now convert to healpix coordinates
                                            healpix_index = hp.ang2pix(nside,θ, ɸ_rot)
                                                    
                                            # Append channel power to ref healpix map
                                            for i in range(len(healpix_index)):
                                                ref_map[healpix_index[i]].append(channel_power[i])
                                            #[ref_map[healpix_index[i]].append(channel_power[i]) for i in range(len(healpix_index))]
                                             
                                            
                                            # Increment pix ounter to keep track of passes in each pix 
                                            for i in healpix_index:
                                                ref_counter[i] += 1
        
            except Exception:
                # Exception message is forwarded from ../decode_rf_data/rf_data.py
                continue

    # Save map arrays to npz file
    np.savez_compressed(f'{out_dir}/{ref}_map_healpix.npz',
            ref_map = ref_map,
            ref_counter = ref_counter
            )


        
if __name__=='__main__':

    parser = argparse.ArgumentParser(description="""
        Project a sat pass onto a healpix map using ephemeris data
        """)
    
    parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
    parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/null_test/',help='Output directory. Default=./../../outputs/null_test/')
    parser.add_argument('--plt_dir', metavar='\b', default='./../../outputs/null_test/pass_plots/',help='Output directory. Default=./../../outputs/null_test/pass_plots/')
    parser.add_argument('--ali_dir', metavar='\b', default='./../../outputs/align_data/',help='Output directory. Default=./../../outputs/align_data/')
    parser.add_argument('--map_dir', metavar='\b', default='../../outputs/sat_channels/window_maps/',help='Satellite channel map. Default=../../outputs/sat_channels/window_maps/')
    parser.add_argument('--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
    parser.add_argument('--noi_thresh', metavar='\b', default=3,help='Noise Threshold: Multiples of MAD. Default=3.')
    parser.add_argument('--sat_thresh', metavar='\b', default=1,help='1 σ threshold to detect sats Default=1.')
    parser.add_argument('--nside', metavar='\b', default=32,help='Healpix Nside. Default = 32')
    parser.add_argument('--plots', metavar='\b', default=False,help='If True, create a gazzillion plots for each sat pass. Default = False')
    
    args = parser.parse_args()
    
    chrono_dir =        args.chrono_dir
    start_date =        args.start_date
    stop_date =         args.stop_date
    out_dir =           args.out_dir
    plt_dir =           args.plt_dir
    ali_dir =           args.ali_dir
    map_dir =           args.map_dir
    noi_thresh =        args.noi_thresh
    sat_thresh =        args.sat_thresh
    nside =             args.nside
    plots =             args.plots

    ref_names=['rf0XX', 'rf0YY', 'rf1XX', 'rf1YY']
    
    if plots == 'True':
        Path(plt_dir).mkdir(parents=True, exist_ok=True)
    
    # Save logs 
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f'{out_dir}/logs_{start_date}_{stop_date}.txt', 'a')

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(proj_ref_healpix, ref_names)

    

