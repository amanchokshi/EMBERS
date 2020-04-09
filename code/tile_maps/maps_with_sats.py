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

map_data = np.load('../../outputs/tile_maps/S08XX_rf0XX_healpix_map.npz', allow_pickle=True)
map_data = {key:map_data[key].item() for key in map_data}

nside = 32

pointings = ['0', '2', '4']

good_sats = [
        25984, 40087, 40091, 
        41179, 41180, 41182,
        41183, 41184, 41185,
        41187, 41188, 41189, 
        28654, 25338, 44387
        ]

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
#    plot_healpix(data_map=ratio_sat_scaled, sub=(1,1,1), cmap=jade, vmin=-30, vmax=0)
#    plt.savefig(f'maps/{s}.png',bbox_inches='tight')
#    #plt.show()
#    plt.close()


#good_map = np.asarray([[] for pixel in range(hp.nside2npix(nside))])
#for sats in good_sats:
#    np.concatenate()

print(len(ratio_sat_data[good_sats[0]]))
print(len(ratio_sat_data[good_sats[1]]))

good_map = np.concatenate((ratio_sat_data[good_sats[0]], ratio_sat_data[good_sats[1]]),axis=1)

good_map_med = [(np.median(i) if i != [] else np.nan ) for i in good_map]
good_map_scaled = np.asarray([(i - np.nanmax(good_map_med[:5000])) for i in good_map_med])
plot_healpix(data_map=good_map_scaled, sub=(1,1,1), cmap=jade, vmin=-30, vmax=0)
plt.savefig(f'maps/good_map.png',bbox_inches='tight')
#plt.show()
plt.close()

#print(ratio_sat_data[41188])

