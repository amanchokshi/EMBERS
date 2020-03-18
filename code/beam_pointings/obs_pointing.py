import sys
import math
import json
import pytz
import json
import argparse
import numpy as np
from pathlib import Path


from pathlib import Path
from astropy.time import Time
from datetime import datetime, timedelta


parser = argparse.ArgumentParser(description="""
        Classify the pointings of 30 min observations
        """)

parser.add_argument('--start_date', metavar='\b', help='Date from which to determine sat ephemeris. Ex: 2019-10-10')
parser.add_argument('--stop_date', metavar='\b', help='Date until which to determine sat ephemeris. Ex: 2019-10-11')
parser.add_argument('--time_zone', metavar='\b', default='Australia/Perth', help='Time zone where data was recorded. Default=Australia/Perth')
parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/beam_pointings/', help='Directory where json metadata files live. Default=./../../outputs/beam_pointings/')
parser.add_argument('--f_name', metavar='\b', default='ultimate_pointing_times.json', help='File name of json to be plotted. Default=ultimate_pointing_times.json')

args = parser.parse_args()
start_date  = args.start_date
stop_date   = args.stop_date
time_zone   = args.time_zone
out_dir     = Path(args.out_dir)
f_name      = Path(args.f_name)



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

# Start and end of 30 min obs in gps time
# Round to nearest int
obs_gps     = np.rint(Time(obs_unix, format='unix').gps)
obs_gps_end = np.rint(Time(obs_unix_end, format='unix').gps)


with open(f'{out_dir}/{f_name}') as table:
    data = json.load(table)
    pointings   = data['grid_pt'] 
    start_gps   = data['start_gps']
    stop_gps    = data['stop_gps'] 
    obs_length  = data['obs_length']

#for i in range(len(start_gps)):
#    print(f'{pointings[i]}: {start_gps[i]}: {stop_gps[i]}: {obs_length[i]}')
#
#print(start_gps[0], start_gps[-1])
#print(obs_gps[0], obs_gps[-1])


def time_filter(s_rise, s_set, times):
    
    # I. sat rises before times, sets within times window
    if (s_rise < times[0] and s_set > times[0] and s_set <= times[-1]):
        occu = (s_set - times[0])/1800


    
    # II. sat rises and sets within times
    elif (s_rise >= times[0] and s_set <= times[-1]):
        occu = (s_set - s_rise)/1800
    
    # III. sat rises within times, and sets after
    elif (s_rise >= times[0] and s_rise < times[-1] and s_set > times[-1]):
        occu = (times[-1] - s_rise)/1800
    
    # IV. sat rises before times and sets after
    elif (s_rise < times[0] and s_set > times[-1]):
        occu = (times[-1] - times[0])/1800
   
    # V. sat completely out of times. Could be on either side
    else:
        occu = None
    
    # intvl = interval
    return occu

point_0 = []
point_2 = []
point_4 = []

for i in range(len(obs_time)):
    obs_window = [obs_gps[i], obs_gps_end[i]]

    for j in range(len(start_gps)):
        occu = time_filter(start_gps[j], stop_gps[j], obs_window)

        if (occu != None and occu >= 0.7):
            if pointings[j] == 0:
                point_0.append(obs_time[i])
            elif pointings[j] == 2:
                point_2.append(obs_time[i])
            else:
                point_4.append(obs_time[i])


# Create dictionary to be saved to json
obs_pointings = {}
obs_pointings['point_0'] = point_0
obs_pointings['point_2'] = point_2
obs_pointings['point_4'] = point_4


with open(f'{out_dir}/obs_pointings.json', 'w') as outfile:
    json.dump(obs_pointings, outfile, indent=4)












