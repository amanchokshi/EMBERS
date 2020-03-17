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
            '--chrono_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json',
            help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')

    parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
    parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
    parser.add_argument('--noi_thresh', metavar='\b', default=3,help='Noise Threshold: Multiples of MAD. Default=3.')
    parser.add_argument('--sat_thresh', metavar='\b', default=3,help='3 σ threshold to detect sats Default=3.')
    parser.add_argument('--nside', metavar='\b', default=32,help='Healpix Nside. Default = 32')
    
    args = parser.parse_args()
    
    start_date        = args.start_date
    stop_date         = args.stop_date
    chan_map          = args.chan_map
    noi_thresh        = args.noi_thresh
    sat_thresh        = args.sat_thresh
    nside             = args.nside
    
    align_dir         = Path(args.align_dir)
    chrono_dir        = Path(args.chrono_dir)
    out_dir           = Path(args.out_dir)

    # Import list of Norad catalogue IDs
    sat_list = [id for id in sat_ids.norad_ids.values()]

    # Tile names
    tile_n  = tile_names()
    refs    = tile_n[:4]
    tiles   = tile_n[4:]

    refx    = refs[::2]
    refy    = refs[1::2]
    tilesx  = tiles[::2]
    tilesy  = tiles[1::2]

    # Read channel map file
    with open(chan_map) as map:
        channel_map = json.load(map)

    # dates: list of days
    # date_time = list of 30 min observation windows
    dates, date_time = time_tree(start_date, stop_date)
    for day in range(len(dates)):
        for window in range(len(date_time[day])):
            for x in tilesx:
                if Path(f'{align_dir}/{dates[day]}/{date_time[day][window]}/rf0XX_{x}_{date_time[day][window]}_aligned.npz').is_file():
                    print(f'rf0XX_{x}_{date_time[day][window]}_aligned.npz exists')
                else:
                    print(f'rf0XX_{x}_{date_time[day][window]}_aligned.npz missing')
            #print(date_time[day][window])
