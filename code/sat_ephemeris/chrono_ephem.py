import os
import math
import json
import pytz
import argparse
import numpy as np

from scipy import interpolate
from astropy.time import Time
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
parser.add_argument('--interp_type', metavar='\b', default='cubic', help='Type of interpolation. Ex: Cubic,Linear. Default=Cubic')
parser.add_argument('--interp_freq', metavar='\b', default=2, help='Frequency at which to interpolate, in Hertz. Must be the same as used in align_data.py. Default=2')
parser.add_argument('--start_date', metavar='\b', help='Date from which to determine sat ephemeris. Ex: 2019-10-10')
parser.add_argument('--stop_date', metavar='\b', help='Date until which to determine sat ephemeris. Ex: 2019-10-11')
parser.add_argument('--time_zone', metavar='\b', default='Australia/Perth', help='Time zone where data was recorded. Default=Australia/Perth')

args = parser.parse_args()

json_dir    = args.json_dir
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


for file in os.listdir(json_dir):
    if file.endswith('.json'):
        f_path = os.path.join(json_dir, file)
        
        with open(f_path) as ephem:
            sat_ephem = json.load(ephem)

            t_array = sat_ephem['time_array']
            # Get start times of all passes
            t_rise = [i[0] for i in t_array]
            unix_rise = Time(t_rise, format='gps').unix
            #print(unix_rise)
    
            
            for t in range(len(obs_time)):
                for rise in unix_rise:
                    if rise > obs_unix[t] and rise < obs_plus[t]:
                        print(obs_time[t], rise)


            s_alt   = sat_ephem['sat_alt']
            s_az    = sat_ephem['sat_az']
            
            #time_array = []
            #sat_alt = []
            #sat_az = []
            #sat_id = []
           
            
            for p in range(len(t_array)):       # p == pass

                # Don't consider passes with less than 4 samples (~ 1 min)
                if len(t_array[p]) > 3:
                   
                    # Convert gps time from Skyfield to unix time
                    # Because rf explorers record data with unix timestamp
                    times_unix = Time(t_array[p], format='gps').unix
                    
                    # Create interpolation functions
                    alt_interp = interpolate.interp1d(times_unix, s_alt[p], kind=interp_type)
                    az_interp  = interpolate.interp1d(times_unix, s_az[p], kind=interp_type)
                    
                    # Determine a integer interval between which to interpolate
                    t_start = math.ceil(times_unix[0])
                    t_stop  = math.floor(times_unix[-1])

                    # Times array, at which to determine alt/az of sat
                    time_interp = list(np.double(np.arange(t_start, t_stop, (1/interp_freq))))
                    #time_array.extend(time_interp)

                    n_steps = len(time_interp)
                    #sat_num = list(np.repeat(int(sat_ephem['sat_id'][0]), n_steps))
                    sat_num = sat_ephem['sat_id'] * n_steps
                    #sat_id.extend(sat_num)
                    
                    #sat_alt.extend(list(alt_interp(time_interp)))
                    #sat_alt.extend(list(alt_interp(time_interp)))
                    sat_az = az_interp(time_interp)
                    sat_az = az_interp(time_interp)

                    #print(min(time_interp))
                    #print(max(time_interp))

                else:
                    pass
                #break
        break              

#for file in os.listdir(json_dir):
#    if file.endswith('.json'):
#        f_path = os.path.join(json_dir, file)
#        
#        with open(f_path) as ephem:
#            sat_ephem = json.load(ephem)
#
#            t_array = sat_ephem['time_array']
#            s_alt   = sat_ephem['sat_alt']
#            s_az    = sat_ephem['sat_az']
#            
#            #time_array = []
#            #sat_alt = []
#            #sat_az = []
#            #sat_id = []
#           
#            
#            for p in range(len(t_array)):       # p == pass
#
#                # Don't consider passes with less than 4 samples (~ 1 min)
#                if len(t_array[p]) > 3:
#                   
#                    # Convert gps time from Skyfield to unix time
#                    # Because rf explorers record data with unix timestamp
#                    times_unix = Time(t_array[p], format='gps').unix
#                    
#                    # Create interpolation functions
#                    alt_interp = interpolate.interp1d(times_unix, s_alt[p], kind=interp_type)
#                    az_interp  = interpolate.interp1d(times_unix, s_az[p], kind=interp_type)
#                    
#                    # Determine a integer interval between which to interpolate
#                    t_start = math.ceil(times_unix[0])
#                    t_stop  = math.floor(times_unix[-1])
#
#                    # Times array, at which to determine alt/az of sat
#                    time_interp = list(np.double(np.arange(t_start, t_stop, (1/interp_freq))))
#                    #time_array.extend(time_interp)
#
#                    n_steps = len(time_interp)
#                    #sat_num = list(np.repeat(int(sat_ephem['sat_id'][0]), n_steps))
#                    sat_num = sat_ephem['sat_id'] * n_steps
#                    #sat_id.extend(sat_num)
#                    
#                    #sat_alt.extend(list(alt_interp(time_interp)))
#                    #sat_alt.extend(list(alt_interp(time_interp)))
#                    sat_az = az_interp(time_interp)
#                    sat_az = az_interp(time_interp)
#
#                    print(min(time_interp))
#                    print(max(time_interp))
#
#                else:
#                    pass
#                break
#        break              

#mask = list(np.argwhere( (np.asarray(time_array)) > obs_unix[0]) & (np.asarray(time_array) < obs_unix[1]) )
#print(list(compress(time_array, mask)))
#print(len(time_array))


## Save collated sat ephem data to one big json file
#with open(f'{json_dir}/ultimate_ephem_list.json', 'w') as outfile:
#    json.dump(chrono_ephem, outfile, indent=4)


