import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad
import sys
sys.path.append('../sat_ephemeris')
from sat_ids import norad_ids
from plot_healpix import plot_healpix 

sys.path.append('../decode_rf_data')
from colormap import spectral, jade, kelp
jade, _ = jade()

map_data = np.load('../../outputs/tile_maps/S08XX_rf0XX_healpix_map.npz', allow_pickle=True)
map_data = {key:map_data[key].item() for key in map_data}

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

ratio_sat_data = {s:[(np.asarray(ratio_map[i])[np.where(np.asarray(sat_map[i]) == s)]).tolist() for i in range(len(ratio_map))] for s in sat_ids}

fig = plt.figure(figsize=(8,10))
#fig.suptitle(f'Healpix Map: {tile}/{ref} @ {p}', fontsize=16)
ratio_sat_med = [(np.median(i) if i != [] else np.nan ) for i in ratio_sat_data[41188]]
ratio_sat_scaled = np.asarray([(i - np.nanmax(ratio_sat_med[:5000])) for i in ratio_sat_med])
plot_healpix(data_map=ratio_sat_scaled, sub=(1,1,1), cmap=jade, vmin=-30, vmax=0)
#plt.savefig(f'{out_dir}/tile_maps/{tile}_{ref}_{p}_map.png',bbox_inches='tight')
plt.show()
plt.close()

#print(ratio_sat_data[41188])

