import sys
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad

sys.path.append('../sat_ephemeris')
from sat_ids import norad_ids
from plot_healpix import plot_healpix 

sys.path.append('../decode_rf_data')
from colormap import jade
jade, _ = jade()


def sort_sat_map(npz_map):
    '''
    Sorts map data by sat id
    '''
    
    f_name, _ = npz_map.name.split('.')
    tile, ref, _, _ = f_name.split('_')

    # pointings at which to make maps
    pointings = ['0', '2', '4']
   
    # read npz map data
    map_data = np.load(npz_map, allow_pickle=True)
    map_data = {key:map_data[key].item() for key in map_data}
   
    # list of all possible satellites
    sat_ids = list(norad_ids.values())
    
    # create a dictionary, with the keys being sat_ids and the values being healpix maps of data from those sats
    # Fist level keys are pointings with dictionaty values
    # Second level keys are satellite ids with healpix map lists
    ratio_sat_data = {}
    ref_sat_data = {}
    tile_sat_data = {}
    time_sat_data = {}

    # loop over pointings for all maps
    for p in pointings:
        ratio_map   = np.asarray(map_data['healpix_maps'][p])
        tile_map    = np.asarray(map_data['tile_maps'][p])
        ref_map     = np.asarray(map_data['ref_maps'][p])
        sat_map     = np.asarray(map_data['sat_map'][p])
        time_map    = np.asarray(map_data['times'][p])

        ratio_sat_data[p]   = {}
        ref_sat_data[p]     = {}
        tile_sat_data[p]    = {}
        time_sat_data[p]    = {}
        
        # loop over all sats
        for s in sat_ids:

            ratio_sat_data[p][s]    = [] 
            ref_sat_data[p][s]      = [] 
            tile_sat_data[p][s]     = [] 
            time_sat_data[p][s]     = [] 
            

            # loop over every healpix pixel
            for i in range(len(ratio_map)):

                # Find subset of data for each sat
                sat_idx = np.where(np.asarray(sat_map[i]) == s)
                ratio_sat_data[p][s].append((np.asarray(ratio_map[i])[sat_idx]).tolist())
                ref_sat_data[p][s].append((np.asarray(ref_map[i])[sat_idx]).tolist())
                tile_sat_data[p][s].append((np.asarray(tile_map[i])[sat_idx]).tolist())
                time_sat_data[p][s].append((np.asarray(time_map[i])[sat_idx]).tolist())

    # Save map arrays to npz file
    tile_data = {'ratio_map':ratio_sat_data, 'ref_map':ref_sat_data, 'tile_map':tile_sat_data, 'time_map':time_sat_data}
    np.savez_compressed(f'{out_dir}/{tile}_{ref}_sat_maps.npz', **tile_data)
    

if __name__=='__main__':
    
    import argparse
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pathlib import Path
    import concurrent.futures
    
    parser = argparse.ArgumentParser(description="""
        Sort the projected healipix data according to satellites.
        This enables us to determine good satellites and also make
        tile maps with a subset of sat data.
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/tile_maps/tile_maps_sats/',help='Output directory. Default=./../../outputs/tile_maps/tile_maps_sats/')
    parser.add_argument('--map_dir', metavar='\b', default='./../../outputs/tile_maps/',help='Output directory. Default=./../../outputs/tile_maps/')
    parser.add_argument('--nside', metavar='\b', type=int,  default=32,help='Healpix Nside. Default = 32')
    
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    map_dir = Path(args.map_dir)
    nside   = args.nside
    
    map_files = [item for item in map_dir.glob('*.npz')]

    Path(f'{out_dir}').mkdir(parents=True, exist_ok=True)
   
   

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(sort_sat_map, map_files)


