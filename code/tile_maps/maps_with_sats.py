import sys
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad

sys.path.append('../sat_ephemeris')
from sat_ids import norad_ids
from plot_healpix import plot_healpix 

sys.path.append('../decode_rf_data')
from colormap import spectral, jade, kelp
jade, _ = jade()

map_data = np.load('../../outputs/tile_maps/10dB_1sigma/S07XX_rf0XX_healpix_map.npz', allow_pickle=True)
map_data = {key:map_data[key].item() for key in map_data}


sat_ids = list(norad_ids.values())

ratio_map   = np.asarray(map_data['healpix_maps'][pointings[0]])
ref_map     = np.asarray(map_data['ref_maps'][pointings[0]])
tile_map    = np.asarray(map_data['tile_maps'][pointings[0]])
sat_map     = np.asarray(map_data['sat_map'][pointings[0]])
time_map    = np.asarray(map_data['times'][pointings[0]])

# create a dictionary, with the keys being sat_ids and the values being healpix maps of data from those sats
ratio_sat_data = {s:[(np.asarray(ratio_map[i])[np.where(np.asarray(sat_map[i]) == s)]).tolist() for i in range(len(ratio_map))] for s in sat_ids}


#for s in sat_ids:
#
#    fig = plt.figure(figsize=(8,10))
#    #fig.suptitle(f'Healpix Map: {tile}/{ref} @ {p}', fontsize=16)
#    ratio_sat_med = [(np.median(i) if i != [] else np.nan ) for i in ratio_sat_data[s]]
#    ratio_sat_scaled = np.asarray([(i - np.nanmax(ratio_sat_med[:5000])) for i in ratio_sat_med])
#    plot_healpix(data_map=ratio_sat_scaled, sub=(1,1,1), cmap=jade)
#    plt.savefig(f'maps_II/{s}.png',bbox_inches='tight')
#    #plt.show()
#    plt.close()


#good_map = [[] for pixel in range(hp.nside2npix(nside))]
#
#for sat in good_sats_II:
#    for p in range(hp.nside2npix(nside)):
#        good_map[p].extend(ratio_sat_data[sat][p])
#
#good_map_med = [(np.median(i) if i != [] else np.nan ) for i in good_map]
#good_map_scaled = np.asarray([(i - np.nanmax(good_map_med[:5000])) for i in good_map_med])
#plot_healpix(data_map=good_map_scaled, sub=(1,1,1), cmap=jade, vmin=np.nanmin(good_map_scaled), vmax=0)
##plt.savefig(f'maps/good_map.png',bbox_inches='tight')
#plt.show()
#plt.close()


if __name__=='__main__':
    
    import argparse
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pathlib import Path
    import concurrent.futures
    
    parser = argparse.ArgumentParser(description="""
        Plot healpix map of reference data
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
   
    # pointings at which to make maps
    pointings = ['0', '2', '4']
   
    # Good sats from which to make plots
    good_sats = [
            25338, 25984, 25985,
            28654, 40086, 40087,
            40091, 41179, 41180,
            41182, 41183, 41184,
            41185, 41187, 41188,
            41189, 44387
            ]

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(map_plots, map_files)

