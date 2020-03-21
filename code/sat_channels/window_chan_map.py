import sys
import json
import argparse
import numpy as np
import healpy as hp
from pathlib import Path
import concurrent.futures
import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad
from sat_channels import time_tree, time_filter,center_of_gravity
from channels_plt import plt_waterfall_pass, plt_channel, sat_plot

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


def power_ephem(
        ref_file,
        chrono_file,
        sat_id,
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
            
            # length of sat pass. Only consider passes longer than 2 minutes
            window_len = w_stop - w_start + 1
            
            if window_len >= 120:
        
                # Slice [crop] the power/times arrays to the times of sat pass
                power_c = power[w_start:w_stop+1, :]
                times_c = times[w_start:w_stop+1]

                # Now we need to find the most suitable channel in this cropped power array
                ##################################

                possible_chans = []
                occu_list = []

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


                                occu_list.append(window_occupancy)
                                
                                if plots == 'True':
                                    # Plots the channel with satellite pass
                                    plt_channel(
                                            f'{plt_dir}/{timestamp}', times_c, channel_power,
                                            s_chan, min_s, max_s, noise_threshold,
                                            arb_thresh,center, cog, cog_thresh,
                                            sat_id, timestamp)
                                
                                possible_chans.append(s_chan)

                    else:
                        continue

                # If channels are identified in the 30 min obs
                n_chans = len(possible_chans)
                
                if n_chans > 0:
                    
                    if plots == 'True':
                        # plot waterfall with sat window and all selected channels highlighted
                        plt_waterfall_pass(
                                f'{plt_dir}/{timestamp}', power, sat_id,
                                w_start, w_stop, possible_chans,
                                timestamp, cmap)
                    
                    # The most lightly channel is one with the highest occupation
                    good_chan = possible_chans[occu_list.index(max(occu_list))]
                    
                    return good_chan

                    channel_power = power_c[:, good_chan]

                    times_sat = np.asarray(norad_ephem["time_array"])
                    
                    alt = np.asarray(norad_ephem["sat_alt"])
                    az  = np.asarray(norad_ephem["sat_az"])
                            
                    if plots == 'True': 
                        # Plot ephemeris of lightly sats present in ephem window of norad_id
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
                                
                                if len(chrono_ephem[s]["sat_alt"][e_0:e_1+1]) >= occ_thresh*len(times_c):
                                    
                                    if chrono_ephem[s]["sat_id"][0] == sat_id:
                                    
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
                        sat_plot(f'{plt_dir}/{timestamp}', plt_ids, sat_id, plt_alt, plt_az, len(plt_ids), timestamp, 'passes')
                    

                    else:
                        return 0

                else:
                    return 0

            else:
                return 0
        
        else:
            return 0
                

def window_chan_map(obs_stamp):

    date, timestamp = obs_stamp

    channel_map = {}

    ref_file = f'{ali_dir}/{date}/{timestamp}/rf0XX_S06XX_{timestamp}_aligned.npz'
    chrono_file = f'{chrono_dir}/{timestamp}.json'
    
    try:
        Path(ref_file).is_file()

        with open(chrono_file) as chrono:
            chrono_ephem = json.load(chrono)

            if chrono_ephem != []:
        
                norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(len(chrono_ephem))]

                if norad_list != []:
                
                    for sat in norad_list:

                        sat_data = power_ephem(
                                ref_file,
                                chrono_file,
                                sat,
                                timestamp
                                )
                        
                        if sat_data != 0:
                            sat_chan = sat_data

                            channel_map[f'{sat}'] = sat_chan

    
    except Exception:
        # Exception message is forwarded from ../decode_rf_data/rf_data.py
        continue

    # Save channel map
    with open(f'{out_dir}/{ref_tile}/channel_data/{int(norad_id)}.json','w') as f: 
        json.dump(sat_data, f, indent=4) 


        
if __name__=='__main__':

    parser = argparse.ArgumentParser(description="""
        Project a sat pass onto a healpix map using ephemeris data
        """)
    
    parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
    parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_channels/',help='Output directory. Default=./../../outputs/sat_channels/')
    parser.add_argument('--plt_dir', metavar='\b', default='./../../outputs/sat_channels/pass_plots',help='Output directory. Default=./../../outputs/sat_channels/pass_plots/')
    parser.add_argument('--ali_dir', metavar='\b', default='./../../outputs/align_data/',help='Output directory. Default=./../../outputs/align_data/')
    parser.add_argument('--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
    parser.add_argument('--noi_thresh', metavar='\b', default=3,help='Noise Threshold: Multiples of MAD. Default=3.')
    parser.add_argument('--sat_thresh', metavar='\b', default=1,help='1 σ threshold to detect sats Default=1.')
    parser.add_argument('--arb_thresh', metavar='\b', default=12,help='Arbitrary Threshold to detect sats. Default=12 dB.')
    parser.add_argument('--occ_thresh', metavar='\b', default=0.80,help='Occupation Threshold of sat in window. Default=0.70')
    parser.add_argument('--cog_thresh', metavar='\b', default=0.05,help='Center of Gravity Threshold to detect sats. Default=0.05')
    parser.add_argument('--plots', metavar='\b', default=False,help='If True, create a gazzillion plots for each sat pass. Default = False')
    
    args = parser.parse_args()
    
    chrono_dir =        args.chrono_dir
    start_date =        args.start_date
    stop_date =         args.stop_date
    out_dir =           args.out_dir
    plt_dir =           args.plt_dir
    ali_dir =           args.ali_dir
    chan_map =          args.chan_map
    noi_thresh =        args.noi_thresh
    sat_thresh =        args.sat_thresh
    arb_thresh =        args.arb_thresh
    occ_thresh =        args.occ_thresh
    cog_thresh =        args.cog_thresh
    nside =             args.nside
    plots =             args.plots

    ref_names=['rf0XX', 'rf0YY', 'rf1XX', 'rf1YY']
    
    if plots == 'True':
        Path(plt_dir).mkdir(parents=True, exist_ok=True)
    
    # Save logs 
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f'{out_dir}/logs_{start_date}_{stop_date}.txt', 'a')

    # Help traverse all 30 min obs b/w start & stop
    dates, date_time = time_tree(start_date, stop_date)
    obs_list = [[dates[d], date_time[d][dt]] for d in range(len(dates)) for dt in range(len(date_time[d]))]

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(proj_ref_healpix, ref_names)

    

