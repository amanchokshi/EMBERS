import os
import sys
import json
import argparse
import sat_ephemeris as se
import sat_ids

import numpy as np
import skyfield as sf
import concurrent.futures
from astropy.time import Time
from skyfield.api import Topos, load


# Skyfield Timescale
ts = load.timescale(builtin=True)

# Position of MWA site in Lat/Lon/Elevation
MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)

parser = argparse.ArgumentParser(description="""
        Code which converts the TLE files downloaded with download_TLE.py
        into satellite ephemeris data: rise time, set time, alt/az arrays
        at a given time cadence. This is saved to a json file which will 
        be used to plot the satellite passes.
        """)
        
parser.add_argument('--sat', metavar='\b', help='The Norad cat ID of satellite. Example: 21576')
parser.add_argument('--tle_dir', metavar='\b', default='./../../outputs/sat_ephemeris/TLE', help='Directory where TLE files are saved. Default=./../../outputs/sat_ephemeris/TLE')
parser.add_argument('--cadence', metavar='\b', default=20, help='Rate at which sat alt/az is computed. Expensive! Default=20s')
parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_ephemeris/ephem_json/', help='Path to output directory. Default=./../../outputs/sat_ephemeris/ephem_json/')


args = parser.parse_args()
sat_name = args.sat
tle_dir = args.tle_dir
cadence = int(args.cadence)
out_dir = args.out_dir

tle_path = f'{tle_dir}/{sat_name}.txt'

sat_list = list(sat_ids.norad_ids.values())


# Save logs 
Path(out_dir).mkdir(parents=True, exist_ok=True)
sys.stdout = open(f'{out_dir}/sat_ephem_logs.txt', 'a')

os.makedirs(os.path.dirname(out_dir), exist_ok=True)

def sat_json(sat_id):
    try:
        sat_ephem = {}
        sat_ephem['sat_id'] = [sat_name]
        sat_ephem['t_rise'] = []
        sat_ephem['t_set'] = []
        sat_ephem['sat_alt'] = []
        sat_ephem['sat_az'] = []
        
        sats, epochs = se.load_tle(tle_path)
        epoch_range = se.epoch_ranges(epochs)
        
        for i in range(len(epoch_range) - 1):   
            t_arr, index_epoch = se.epoch_time_array(epoch_range, i, cadence)
            passes, alt, az = se.sat_pass(sats, t_arr, index_epoch) 
            print(passes)
            for pass_index in passes:
                t_rise, t_set, sat_alt, sat_az = se.ephem_data(t_arr, pass_index, alt, az)
        
                sat_ephem['t_rise'].append(t_rise)
                sat_ephem['t_set'].append(t_set)
                sat_ephem['sat_alt'].append(sat_alt)
                sat_ephem['sat_az'].append(sat_az)
        
        with open(f'{out_dir}/{sat_name}.json', 'w') as outfile:
            json.dump(sat_ephem, outfile) 
        
        return f'Saved {sat_name}.json'
    except Exception:
        return f'ERROR! Couldn\'t save {sat_name}.json.'

with concurrent.futures.ProcessPoolExecutor() as executor:
    results = executor.map(sat_json, sat_list)

for result in results:
    print(result)

