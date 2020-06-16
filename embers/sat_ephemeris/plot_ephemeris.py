import os
import json
import sat_ids
import concurrent.futures
from itertools import repeat
from sat_ephemeris import sat_plot


def plot_ephem(json_dir, out_dir, sat_id):
    try:
        with open(f'{json_dir}/{sat_id}.json', 'r') as ephem:
            sat_ephem = json.load(ephem)
            altitude = sat_ephem['sat_alt']
            azimuth = sat_ephem['sat_az']
            num_passes = len(sat_ephem['time_array'])
            
            plt = sat_plot(sat_id, altitude, azimuth, num_passes)
            plt.savefig(f'{out_dir}/{sat_id}.png')
            plt.close()

    except Exception:
        print(f'Could not plot {sat_id}.png. Source json file missing ')

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description=""" 
            Plot Satellite passes from json ephemeris files.
            """)
    
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_ephemeris/ephem_plots/', help='Path to output directory. Default=./../../outputs/sat_ephemeris/ephem_plots/')
    parser.add_argument('--json_dir', metavar='\b', default='./../../outputs/sat_ephemeris/ephem_json/', help='Path to ephem_json directory. Default=./../../outputs/sat_ephemeris/ephem_json/')
    
    args = parser.parse_args()
    
    out_dir =   args.out_dir
    json_dir =  args.json_dir
    
    # Create output directory
    os.makedirs(os.path.dirname(out_dir), exist_ok=True) 
    
    # Load list of satellites (Nodar ids)
    sat_list = list(sat_ids.norad_ids.values())
    
    # Parellization Magic Here!
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(
                plot_ephem, 
                repeat(json_dir),
                repeat(out_dir),
                sat_list)
    
