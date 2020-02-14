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

# Number of days that the date range spans
n_days = (t_stop - t_start).days

# YYYY-MM-DD-HH:MM format
obs_time = []

# Start of half hour obs in unix time
obs_unix = []

# End of half hour obs in unix time
obs_unix_end = []

# The +1 makes the date ranges inclusive
for i in range(n_days+1):
    day = t_start + timedelta(days=i)
    date = day.strftime('%Y-%m-%d')
    
    # loop over 48x30 minute obs in a day
    for j in range(48):
        t_delta = datetime.strptime(date,'%Y-%m-%d') + timedelta(minutes=30*j)
        # convert t_delta to a readable string YYYY-MM-DD-HH:MM
        d_time = t_delta.strftime('%Y-%m-%d-%H:%M')
        
        # Convert from naive local time to utc aware time
        utc_delta = local.localize(t_delta, is_dst=None).astimezone(pytz.utc)
        
        # convert to a unix timestamp, used within rf explorer data files
        utc_unix = utc_delta.timestamp()
        # time at end of half hour window
        utc_unix_end = utc_unix + (30 * 60)
        
        obs_time.append(d_time)
        obs_unix.append(utc_unix)
        obs_unix_end.append(utc_unix_end)


def interp_ephem(t_array, s_alt, s_az, interp_type, interp_freq):
    '''Interpolates satellite ephemeris - time, alt, az
    
    Interpolate sat ephem to same freq as align data. This ensures
    that each point in the data, will have an corresponding ephem
    point.

    Args:
        t_array:        Time array of on satellite pass
        s_alt, s_az:    Satellite alt, az at the t_array
        interp_type:    Type of interpolation. Ex: Cubic,Linear. Default=Cubic'
        interp_freq:    Frequency at which to interpolate, in Hertz. 
        
    Returns:
        time_interp:    Interpolated time array
        sat_alt:        Interpolated alt array
        sat_az:         Interpolated az array
    '''

    # Don't consider passes with less than 4 samples (~ 1 min = 3*20s) 
    if len(t_array) > 3:
        
        # Create interpolation functions. Math functions, not Python!
        alt_interp = interpolate.interp1d(t_array, s_alt, kind=interp_type)
        az_interp  = interpolate.interp1d(t_array, s_az, kind=interp_type)
        
        # Makes start and end times clean integers
        # Also ensures that the interp range is inclusive of data points
        start = math.ceil(t_array[0])
        stop = math.floor(t_array[-1])

        # Create time array, at which to evaluate alt/az of sat
        time_interp = list(np.double(np.arange(start, stop, (1/interp_freq))))
    
        sat_alt = list(alt_interp(time_interp))
        sat_az = list(az_interp(time_interp))

        return(time_interp, sat_alt, sat_az)
    
    else:
        pass







json_path = '../../outputs/sat_ephemeris/ephem_json/21576.json'

with open(json_path) as ephem:
    sat_ephem = json.load(ephem)
    
    # Extract data from json dictionary
    t_array = sat_ephem['time_array']
    s_alt   = sat_ephem['sat_alt']
    s_az    = sat_ephem['sat_az']
    s_id  = sat_ephem['sat_id']

    print(len(t_array))

for pass_idx in range(len(t_array)):
    time_interp, sat_alt, sat_az = interp_ephem(
            t_array[pass_idx],
            s_alt[pass_idx],
            s_az[pass_idx],
            interp_type,
            interp_freq)

    for obs_int in range(len(obs_unix)):

        # Case I: Satpass occurs completely within the 30min observation
        if obs_unix[obs_int] <= time_interp[0] and obs_unix_end[obs_int] >= time_interp[-1]:
            print(f'{pass_idx}: I {obs_time[obs_int]}')

        # Case II: Satpass begins before the obs, but ends within it
        elif obs_unix[obs_int] >= time_interp[0] and obs_unix[obs_int] <= time_interp[-1] and obs_unix_end[obs_int] >= time_interp[-1]:
            print(f'{pass_idx}: II {obs_time[obs_int]}')

        # Case III: Satpass begins within the obs, but ends after it
        elif obs_unix_end[obs_int] >= time_interp[0] and obs_unix_end[obs_int] <= time_interp[-1] and obs_unix[obs_int] <= time_interp[0]:
            print(f'{pass_idx}: III {obs_time[obs_int]}')
    
    #break # breaks loop after first pass in ephem.json
  
    



