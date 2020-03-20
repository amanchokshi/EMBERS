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


def power_ephem(
        ali_file,
        chrono_file,
        sat_id,
        sat_chan,
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

            # Only continue, if sat passes are longer then 120 sec
            if (w_stop - w_start) >= 120:
        
                # Slice [crop] the ref/tile/times arrays to the times of sat pass and extract sat_chan
                ref_c = ref_p[w_start:w_stop+1, sat_chan]
                tile_c = tile_p[w_start:w_stop+1, sat_chan]
                
                alt = np.asarray(norad_ephem["sat_alt"])
                az  = np.asarray(norad_ephem["sat_az"])

                # Apply noise criteria. In the window, where are ref_power and tile power
                # above their respective thresholds?
                if np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0] is not []: 
                    good_ref    = ref_c[np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]]
                    good_tile   = tile_c[np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]]
                    good_alt    = alt[np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]]
                    good_az     = az[np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]]
            
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

    # Initialize an empty dictionary for tile data
    # The map is list of length 12288 of empty lists to append pixel values to
    # The counter is an array of zeros, to increment when a new values is added
    tile_data = {'healpix_maps':{p:[[] for pixel in range(hp.nside2npix(nside))] for p in pointings}, 'healpix_counters':{p:np.zeros(hp.nside2npix(nside)) for p in pointings}}
    
    for day in range(len(dates)):

        for window in range(len(date_time[day])):
            timestamp = date_time[day][window]
            
            # Check if at timestamp, reciever was pointed to 0,2,4 gridpointing
            if ((timestamp in point_0) or (timestamp in point_2) or (timestamp in point_4)):
            
                # pointing at timestamp
                point = check_pointing(timestamp, point_0, point_2, point_4)

                #print(timestamp, point)
               
                ali_file = Path(f'{align_dir}/{dates[day]}/{timestamp}/{ref}_{tile}_{timestamp}_aligned.npz')
                
                # check if file exists
                if ali_file.is_file():

                    #print(f'Exists {ref}_{tile}_{timestamp}_aligned.npz')
                    
                    # Chrono Ephemeris file
                    chrono_file = Path(f'{chrono_dir}/{timestamp}.json')

                    with open(chrono_file) as chrono:
                        chrono_ephem = json.load(chrono)
                    
                        if chrono_ephem != []:
                
                            norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(len(chrono_ephem))]
                            
                            for sat in list(channel_map.keys()):
                                                
                                if int(sat) in norad_list and norad_list != []:

                                    chans = channel_map[sat]
                                    
                                    for chan_num in chans:

                                        sat_data = power_ephem(
                                                ali_file,
                                                chrono_file,
                                                int(sat),
                                                chan_num)

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
                                            #[ref_map[healpix_index[i]].append(channel_power[i]) for i in range(len(healpix_index))]
                                             
                                            
                                            # Increment pix ounter to keep track of passes in each pix 
                                            for i in healpix_index:
                                                tile_data['healpix_counters'][f'{point}'][i] += 1
                            
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
    #sys.stdout = open(f'{out_dir}/logs_{start_date}_{stop_date}.txt', 'a')
    
    
#    for tile_pair in tile_pairs:
#        project_tile_healpix(tile_pair)
#        break
         
    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(project_tile_healpix, tile_pairs)
    


