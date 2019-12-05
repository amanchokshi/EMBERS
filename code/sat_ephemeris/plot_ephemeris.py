import os
import json
import argparse
import sat_ids
from sat_ephemeris import sat_plot


parser = argparse.ArgumentParser(description=""" 
        Plot Satellite passes from json ephemeris files.
        """)

#parser.add_argument('--sat', metavar='\b', help='Norad cat ID of satellite. Ex:23545')
parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_ephemeris/ephem_plots/', help='Path to output directory. Default=./../../outputs/sat_ephemeris/ephem_plots/')
parser.add_argument('--json_dir', metavar='\b', default='./../../outputs/sat_ephemeris/ephem_json/', help='Path to ephem_json directory. Default=./../../outputs/sat_ephemeris/ephem_json/')

args = parser.parse_args()
#sat_id = args.sat
out_dir = args.out_dir
json_dir = args.json_dir

# Create output directory
os.makedirs(os.path.dirname(out_dir), exist_ok=True) 

# Load list of satellites (Nodar ids)
sat_list = list(sat_ids.norad_ids.values())

for sat_id in sat_list:
    try:
        with open(f'{json_dir}/{sat_id}.json', 'r') as ephem:
            sat_ephem = json.load(ephem)
            altitude = sat_ephem['sat_alt']
            azimuth = sat_ephem['sat_az']
            num_passes = len(sat_ephem['t_rise'])
            
            plt = sat_plot(sat_id, altitude, azimuth, num_passes)
            plt.savefig(f'{out_dir}/{sat_id}.png')
            plt.close()

    except Exception:
        print(f'Could not plot {sat_id}.png. Source json file missing ')



