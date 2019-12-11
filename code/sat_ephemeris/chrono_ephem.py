import os
import math
import json
import pytz
import argparse
import numpy as np

from pathlib import Path
from scipy import interpolate
from itertools import compress
from datetime import datetime, timedelta

parser = argparse.ArgumentParser(description="""
        Collates satellite pass data from all ephem json files.
        The native skyfiled gps timestamps are converted to unix
        timestamps to match the output of the rf explorers. The 
        alt, az data is interpolated to match the cadence of 
        align_data.py. The passes are sorted based on rise time 
        and saved to a new file - ultimate_ephem_list.json
        """)

parser.add_argument('--json_dir', metavar='\b', default='./../../outputs/sat_ephemeris/ephem_json/', help='Directory where ephem json files live. Default=./../../outputs/sat_ephemeris/ephem_json/')
parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json/', help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
parser.add_argument('--interp_type', metavar='\b', default='cubic', help='Type of interpolation. Ex: Cubic,Linear. Default=Cubic')
parser.add_argument('--interp_freq', metavar='\b', default=2, help='Frequency at which to interpolate, in Hertz. Must be the same as used in align_data.py. Default=2')
parser.add_argument('--start_date', metavar='\b', help='Date from which to determine sat ephemeris. Ex: 2019-10-10')
parser.add_argument('--stop_date', metavar='\b', help='Date until which to determine sat ephemeris. Ex: 2019-10-11')
parser.add_argument('--time_zone', metavar='\b', default='Australia/Perth', help='Time zone where data was recorded. Default=Australia/Perth')

args = parser.parse_args()

json_dir    = args.json_dir
out_dir    = args.out_dir
interp_type = args.interp_type
interp_freq = args.interp_freq
start_date  = args.start_date
stop_date   = args.stop_date
time_zone   = args.time_zone


# Time stuff
# The time input is in local time. As in Austraila/Perth
local = pytz.timezone(time_zone)

t_start = datetime.strptime(start_date, '%Y-%m-%d')
t_stop = datetime.strptime(stop_date, '%Y-%m-%d')
n_days = (t_stop - t_start).days

obs_time = []
obs_unix = []
obs_plus = []

for i in range(n_days+1):
    day = t_start + timedelta(days=i)
    date = day.strftime('%Y-%m-%d')
    
    for j in range(48):
        t_delta = datetime.strptime(date,'%Y-%m-%d') + timedelta(minutes=30*j)
        d_time = t_delta.strftime('%Y-%m-%d-%H:%M')
        
        # Convert from naive local time to utc aware time
        utc_delta = local.localize(t_delta, is_dst=None).astimezone(pytz.utc)
        
        # convert to a unix timestamp, used within rf explorer data files
        utc_unix = utc_delta.timestamp()
        # time at end of half hour window
        utc_plus = utc_unix + (30 * 60)
        
        obs_time.append(d_time)
        obs_unix.append(utc_unix)
        obs_plus.append(utc_plus)



def interp_ephem(start, stop, pass_idx, t_array, s_alt, s_az, interp_type, interp_freq):
    '''Interpolates satellite ephemeris - time, alt, az - to a given frequency'''

    # Don't consider passes with less than 4 samples (~ 1 min)
    if len(t_array[pass_idx]) > 3:
        
        
        # Create interpolation functions
        alt_interp = interpolate.interp1d(t_array[pass_idx], s_alt[pass_idx], kind=interp_type)
        az_interp  = interpolate.interp1d(t_array[pass_idx], s_az[pass_idx], kind=interp_type)
        
        # Makes start ann end times clean integers
        # Also ensures that the interp range is inclusive of data points
        start = math.ceil(start)
        stop = math.floor(stop)

        # Determine a integer interval between which to interpolate
        # Times array, at which to determine alt/az of sat
        time_interp = list(np.double(np.arange(start, stop, (1/interp_freq))))
    
        sat_alt = list(alt_interp(time_interp))
        sat_az = list(az_interp(time_interp))

        return(time_interp, sat_alt, sat_az)
    
    else:
        pass


# loop through sat ephem json files
for file in os.listdir(json_dir):
    if file.endswith('.json'):
        f_path = os.path.join(json_dir, file)
       
        # Work with one sat at a time
        with open(f_path) as ephem:
            sat_ephem = json.load(ephem)
            
            # Extract data from json dictionary
            t_array = sat_ephem['time_array']
            s_alt   = sat_ephem['sat_alt']
            s_az    = sat_ephem['sat_az']
            sat_id  = sat_ephem['sat_id']
            
            
            # Get rise and set times of all passes
            t_rise = [i[0] for i in t_array]
            t_set = [i[-1] for i in t_array]


            # iterate over all the half hour observations 
            for t in range(len(obs_time)):
                
                # Empty lists to hold interpolated data
                time_array = []
                sat_alt = []
                sat_az = []
             
                # iterate over all passes, and see whether any intersect with particular half hour obs
                for idx in range(len(t_rise)):
                    
                    # Case I - pass starts before obs begins, but ends within the observation
                    if t_rise[idx] <= obs_unix[t] and t_set[idx] <= obs_plus[t] and t_set[idx] > obs_unix[t]:
                        tm_start = obs_unix[t]
                        tm_stop = t_set[idx]

                    # Case II - pass starts and ends within the observation
                    elif t_rise[idx] >= obs_unix[t] and t_set[idx] <= obs_plus[t]:
                        tm_start = t_rise[idx]
                        tm_stop = t_set[idx]
                    
                    # Case III - pass starts within obs, but ends after the observation ends
                    elif t_rise[idx] >= obs_unix[t] and t_rise[idx] < obs_plus[t] and t_set[idx] >= obs_plus[t]:
                        tm_start = t_rise[idx]
                        tm_stop = obs_plus[t]
                   
                        # call the interp_ephem function
                        tm_interp, st_alt, st_az = interp_ephem(
                                tm_start,
                                tm_stop,
                                idx,
                                t_array,
                                s_alt,
                                s_az,
                                interp_type,
                                interp_freq)
                        
                        # append relevant interpolated data to lists
                        time_array.append(tm_interp)    
                        sat_alt.append(st_alt)
                        sat_az.append(st_az)

                    else:
                        pass

                
                s_ephem = {}
                s_ephem['sat_id'] = sat_id
                s_ephem['time_array'] = time_array
                s_ephem['sat_alt'] = sat_alt
                s_ephem['sat_az'] = sat_az

                Path(out_dir).mkdir(parents=True, exist_ok=True)

                if time_array != []:
                    with open(f'{out_dir}/{obs_time[t]}.json', 'w') as outfile:
                        json.dump(s_ephem, outfile, indent=4)

        # TODO figure out how to loop over sats, and save to same file
        break
           
           


