import sys
import json
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
from pathlib import Path
import concurrent.futures
from scipy.stats import median_absolute_deviation as mad

sys.path.append('../sat_ephemeris')
import sat_ids

sys.path.append('../sat_channels')
from sat_channels import time_tree, savgol_interp, time_filter

# Custom spectral colormap
sys.path.append('../decode_rf_data')
from rf_data import tile_names
from colormap import spectral
cmap = spectral()


def check_pointing(timestamp, point_0, point_2, point_4):
    '''Check if timestamp is at pointing 0, 2, 4'''
    if timestamp in point_0:
        point = 0
    elif timestamp in point_2:
        point = 2
    else:
        point = 4
    
    return point

def read_aligned(ali_file=None):
    '''Read aligned data npz file'''
    paired_data = np.load(ali_file, allow_pickle=True)
    
    ref_p   = paired_data['ref_p_aligned']
    tile_p  = paired_data['tile_p_aligned']
    times   = paired_data['time_array']

    return [ref_p, tile_p, times]


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
        out_dir, times, ref_c, tile_c,
        ref_noise, tile_noise, chan_num,
        sat_id, pointing, date, point):

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

    ylims = [min(ref_c), max(ref_c), min(tile_c), max(tile_c)]
    min_s = min(ylims)
    max_s = max(ylims)
   
    fig = plt.figure(figsize=(8,6))
    fig.suptitle(f'Satellite [{sat_id}] Pass @ {date} in Channel: [{chan_num}]', y=0.98)

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(times, ref_c, linestyle='-', linewidth=2, alpha=1.0, color='#729d39', label='ref')
    ax1.fill_between(times, ref_c, color='#729d39', alpha=0.7)
    ax1.axhline(ref_noise,linestyle='-', linewidth=2, color='#36622b', label=f'Ref Cut: {ref_noise:.2f} dBm')
    ax1.axhspan(-1, ref_noise, color='#36622b', alpha=0.7)
    ax1.set_ylim([min(ref_c)-1, max(ref_c)+1])
    ax1.set_ylabel('Power [dBm]')
    ax1.set_xticklabels([])
    
    leg = ax1.legend(frameon=True)
    leg.get_frame().set_facecolor('grey')
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)
    
    
    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(times, tile_c, linestyle='-', linewidth=2, alpha=1.0, color='#ff5656', label='tile')
    ax2.fill_between(times, tile_c, color='#ff5656', alpha=0.7)
    ax2.axhline(tile_noise,linestyle='-', linewidth=2, color='#970747', label=f'Tile Cut: {tile_noise:.2f} dBm')
    ax2.axhspan(-1, tile_noise, color='#970747', alpha=0.7)
    ax2.set_ylim([min(tile_c)-1, max(tile_c)+1])
    ax2.set_ylabel('Power [dBm]')
    ax2.set_xlabel('Time [s]')
    
    leg = ax2.legend(frameon=True)
    leg.get_frame().set_facecolor('grey')
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.savefig(f'{out_dir}/{date}_{sat_id}_{chan_num}_channel.png')
    plt.close()


def power_ephem(
        ref, tile,
        ali_file,
        chrono_file,
        sat_id,
        sat_chan,
        point,
        arb_thresh,
        timestamp
        ):

    '''Create power, alt, az arrays at constant cadence'''

    # Read .npz aligned file
    ref_p, tile_p, times = read_aligned(ali_file=ali_file)

    # Scale noise floor to zero and determine noise threshold
    ref_p, ref_noise = noise_floor(sat_thresh, noi_thresh, ref_p)
    tile_p, tile_noise = noise_floor(sat_thresh, noi_thresh, tile_p)
    
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

            # Slice [crop] the ref/tile/times arrays to the times of sat pass and extract sat_chan
            ref_c = ref_p[w_start:w_stop+1, sat_chan]
            tile_c = tile_p[w_start:w_stop+1, sat_chan]
            times_c = times[w_start:w_stop+1]
            
            alt = np.asarray(norad_ephem["sat_alt"])
            az  = np.asarray(norad_ephem["sat_az"])

            if ((max(ref_c) >= arb_thresh) and (max(tile_c) >= arb_thresh)):

                # Apply noise criteria. In the window, where are ref_power and tile power
                # above their respective thresholds?
                if np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0].size != 0: 
                    good_ref    = ref_c[np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]]
                    good_tile   = tile_c[np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]]
                    good_alt    = alt[np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]]
                    good_az     = az[np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]]

                    if plots == 'True':
                        plt_channel(f'{plt_dir}/{tile}_{ref}/{point}', times_c, ref_c, tile_c, ref_noise, tile_noise, sat_chan, sat_id, point, timestamp, point)
                
                    return [good_ref, good_tile, good_alt, good_az]

                else:
                    return 0
            
            else:
                return 0

        else:
            return 0


def project_tile_healpix(tile_pair):

    ref, tile = tile_pair

    pointings = ['0','2','4']
    
    if plots == 'True':
        Path(f'{plt_dir}/{tile}_{ref}/0').mkdir(parents=True, exist_ok=True)
        Path(f'{plt_dir}/{tile}_{ref}/2').mkdir(parents=True, exist_ok=True)
        Path(f'{plt_dir}/{tile}_{ref}/4').mkdir(parents=True, exist_ok=True)

    # Initialize an empty dictionary for tile data
    # The map is list of length 12288 of empty lists to append pixel values to
    # The counter is an array of zeros, to increment when a new values is added
    # keep track of which satellites contributed which data
    tile_data = {'healpix_maps':{p:[[] for pixel in range(hp.nside2npix(nside))] for p in pointings}, 'sat_map':{p:[[] for pixel in range(hp.nside2npix(nside))] for p in pointings}}
    
    for day in range(len(dates)):

        for window in range(len(date_time[day])):
            timestamp = date_time[day][window]
            
            # Check if at timestamp, reciever was pointed to 0,2,4 gridpointing
            if ((timestamp in point_0) or (timestamp in point_2) or (timestamp in point_4)):
            
                # pointing at timestamp
                point = check_pointing(timestamp, point_0, point_2, point_4)

                ali_file = Path(f'{align_dir}/{dates[day]}/{timestamp}/{ref}_{tile}_{timestamp}_aligned.npz')

                # check if file exists
                if ali_file.is_file():

                    # Chrono and map Ephemeris file
                    chrono_file = Path(f'{chrono_dir}/{timestamp}.json')
                    channel_map = f'{map_dir}/{timestamp}.json'

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
                                                    ref, tile,
                                                    ali_file,
                                                    chrono_file,
                                                    sat,
                                                    chan,
                                                    point,
                                                    arb_thresh,
                                                    timestamp
                                                    )

                                            if sat_data != 0:
                                            
                                                ref_power, tile_power, alt, az = sat_data

                                                # DIVIDE OUT SATS
                                                # Here we divide satellite signal by reference signal
                                                # In log space, this is subtraction
                                                pass_power = tile_power - ref_power

                                                # Altitude is in deg while az is in radians
                                                # convert alt to radians
                                                alt = np.radians(alt)
                                                az  = np.asarray(az)

                                                # To convert from Alt/Az to θ/ϕ spherical coordinates
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
                                                    tile_data['healpix_maps'][f'{point}'][healpix_index[i]].append(pass_power[i])
                                                 
                                                
                                                # Keep track of sats in each healpix pixel
                                                for i in range(len(healpix_index)):
                                                    tile_data['sat_map'][f'{point}'][healpix_index[i]].append(sat)
                            
                else:
                    print(f'Missing {ref}_{tile}_{timestamp}_aligned.npz')
                    continue

    # Save map arrays to npz file
    np.savez_compressed(f'{out_dir}/{tile}_{ref}_healpix_map.npz', **tile_data)

if __name__=='__main__':

    parser = argparse.ArgumentParser(description="""
        Project a sat pass onto a healpix map using ephemeris data
        """)
    
    parser.add_argument(
            '--align_dir', metavar='\b',default='../../outputs/align_data/',
            help='Dir where agligned data date is saved. Default:../../outputs/align_data')

    parser.add_argument(
            '--out_dir', metavar='\b', default='./../../outputs/tile_maps/',
            help='Output directory. Default=./../../outputs/tile_maps/')
    
    parser.add_argument(
            '--plt_dir', metavar='\b', default='./../../outputs/tile_maps/pass_plots/',
            help='Output directory. Default=./../../outputs/tile_maps/pass_plots/')

    parser.add_argument(
            '--chan_map', metavar='\b', default='../../data/channel_map.json',
            help='Satellite channel map. Default=../../data/channel_map.json')
    
    parser.add_argument(
            '--obs_point', metavar='\b', default='../../outputs/beam_pointings/obs_pointings.json',
            help='Observation pointing lists. Default=../../outputs/beam_pointings/obs_pointings.json')

    parser.add_argument(
            '--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',
            help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
    
    parser.add_argument(
            '--map_dir', metavar='\b', default='../../outputs/sat_channels/window_maps/',
            help='Satellite channel map. Default=../../outputs/sat_channels/window_maps/')

    parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
    parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
    parser.add_argument('--noi_thresh', metavar='\b', default=3,help='Noise Threshold: Multiples of MAD. Default=3.')
    parser.add_argument('--sat_thresh', metavar='\b', default=1,help='3 σ threshold to detect sats Default=1.')
    parser.add_argument('--arb_thresh', metavar='\b', default=10,help='Arbitrary Threshold to detect sats. Default=10 dB.')
    parser.add_argument('--nside', metavar='\b', default=32,help='Healpix Nside. Default = 32')
    parser.add_argument('--plots', metavar='\b', default=False,help='If True, create a gazzillion plots for each sat pass. Default = False')
    
    args = parser.parse_args()
    
    start_date      = args.start_date
    stop_date       = args.stop_date
    chan_map        = args.chan_map
    obs_point       = args.obs_point
    noi_thresh      = args.noi_thresh
    sat_thresh      = args.sat_thresh
    arb_thresh      = args.arb_thresh
    nside           = args.nside
    plots           = args.plots
    
    align_dir       = Path(args.align_dir)
    chrono_dir      = Path(args.chrono_dir)
    out_dir         = Path(args.out_dir)
    plt_dir         = Path(args.plt_dir)
    map_dir         = Path(args.map_dir)

    # Import list of Norad catalogue IDs
    sat_list = [id for id in sat_ids.norad_ids.values()]

    # Tile names
    refs    = tile_names()[:4]
    tiles   = tile_names()[4:]
    
    # All relevant tile pairs
    tile_pairs = []
    for ref in refs:
        if 'XX' in ref:
            for tile in [t for t in tiles if 'XX' in t]:
                tile_pairs.append([ref, tile])
        else:
            for tile in [t for t in tiles if 'YY' in t]:
                tile_pairs.append([ref,tile])


    # Read channel map file
    with open(chan_map) as map:
        channel_map = json.load(map)

    # Read observation pointing list
    with open(obs_point) as point:
        obs_p = json.load(point)
        point_0 = obs_p['point_0'] 
        point_2 = obs_p['point_2'] 
        point_4 = obs_p['point_4']




    # dates: list of days
    # date_time = list of 30 min observation windows
    dates, date_time = time_tree(start_date, stop_date)

    # Save logs 
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f'{out_dir}/logs_{start_date}_{stop_date}.txt', 'a')
   
#    for tile_pair in tile_pairs:
#        project_tile_healpix(tile_pair)
#        break
         
    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(project_tile_healpix, tile_pairs)
    


