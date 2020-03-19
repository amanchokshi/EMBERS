import sys
import json
import argparse
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
    noise_threshold = μ_noise + noi_thresh*σ_noise
    
    # scale the data so that it has zero median
    data = data - μ_noise
    
    return (data, noise_threshold)


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
            '--chan_map', metavar='\b', default='../../data/channel_map.json',
            help='Satellite channel map. Default=../../data/channel_map.json')
    
    parser.add_argument(
            '--obs_point', metavar='\b', default='../../outputs/beam_pointings/obs_pointings.json',
            help='Observation pointing lists. Default=../../outputs/beam_pointings/obs_pointings.json')

    parser.add_argument(
            '--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',
            help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')

    parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
    parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
    parser.add_argument('--noi_thresh', metavar='\b', default=3,help='Noise Threshold: Multiples of MAD. Default=3.')
    parser.add_argument('--sat_thresh', metavar='\b', default=1,help='3 σ threshold to detect sats Default=1.')
    parser.add_argument('--nside', metavar='\b', default=32,help='Healpix Nside. Default = 32')
    
    args = parser.parse_args()
    
    start_date      = args.start_date
    stop_date       = args.stop_date
    chan_map        = args.chan_map
    obs_point       = args.obs_point
    noi_thresh      = args.noi_thresh
    sat_thresh      = args.sat_thresh
    nside           = args.nside
    
    align_dir       = Path(args.align_dir)
    chrono_dir      = Path(args.chrono_dir)
    out_dir         = Path(args.out_dir)

    # Import list of Norad catalogue IDs
    sat_list = [id for id in sat_ids.norad_ids.values()]

    # Tile names
    tile_n  = tile_names()
    refs    = tile_n[:4]
    tiles   = tile_n[4:]
    
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
    
    for day in range(len(dates)):
        
        for window in range(len(date_time[day])):
            timestamp = date_time[day][window]
            point = check_pointing(timestamp, point_0, point_2, point_4)
            
            if ((timestamp in point_0) or (timestamp in point_2) or (timestamp in point_4)):
               
                for pair in tile_pairs:
                    ref, tile = pair
                    f = Path(f'{align_dir}/{dates[day]}/{timestamp}/{ref}_{tile}_{timestamp}_aligned.npz')
                    try:
                        point = check_pointing(timestamp, point_0, point_2, point_4)
                        
                        # Read .npz aligned file
                        ref_p, tile_p, times = read_aligned(ali_file=f)
                       
                        # Scale noise floor to zero and determine noise threshold
                        ref_p, ref_noise = noise_floor(sat_thresh, noi_thresh, ref_p)
                        tile_p, tile_noise = noise_floor(sat_thresh, noi_thresh, tile_p)

                    except FileNotFoundError as e:
                        print(e)




