import os
import json
import argparse
from sat_ephemeris import sat_plot


parser = argparse.ArgumentParser(description=""" 
        Plot Satellite passes from json ephemeris files.
        """)

parser.add_argument('--sat', metavar='\b', help='Norad cat ID of satellite. Ex:23545')

args = parser.parse_args()
sat_id = args.sat

# Create output directory
os.makedirs(os.path.dirname('./outputs/ephem_plots/'), exist_ok=True) 


with open('./outputs/ephem_json/{}.json'.format(sat_id), 'r') as ephem:
    sat_ephem = json.load(ephem)
    altitude = sat_ephem['sat_alt']
    azimuth = sat_ephem['sat_az']
    num_passes = len(sat_ephem['t_rise'])
    
    plt = sat_plot(sat_id, altitude, azimuth, num_passes)
    plt.savefig('')
    plt.show()



