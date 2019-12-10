import os
import math
import json
import pytz
import argparse
import numpy as np

from scipy import interpolate
from astropy.time import Time
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

#time_array  = []
#t_rise      = []
#sat_id      = []
#alt         = []
#az          = []


# Time stuff
# The time input is in local time. As in Austraila/Perth
local = pytz.timezone(time_zone)

t_start = datetime.strptime(start_date, '%Y-%m-%d')
t_stop = datetime.strptime(stop_date, '%Y-%m-%d')
n_days = (t_stop - t_start).days

obs_times = []

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
        
        tm = [d_time, utc_unix, utc_plus]
        obs_times.append(tm)


for obs in obs_times:
    date = obs[0]
    unix_start = obs[1]
    unix_stop = obs[2]

    sat_name = []
    time_unix = []
    alt = []
    az = []



    for file in os.listdir(json_dir):
        if file.endswith('.json'):
            f_path = os.path.join(json_dir, file)
            print(f_path)
            
            with open(f_path) as ephem:
                sat_ephem = json.load(ephem)
    
    
                t_array = sat_ephem['time_array']
                s_alt   = sat_ephem['sat_alt']
                s_az    = sat_ephem['sat_az']
               
                
                for p in range(len(t_array)):
    
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
    
                        n_steps = len(time_interp)
                        #sat_num = list(np.repeat(int(sat_ephem['sat_id'][0]), n_steps))
                        sat_num = sat_ephem['sat_id'] * n_steps
                        
                        sat_alt = list(alt_interp(time_interp))
                        sat_az = list(az_interp(time_interp))
    
    
                        #az.append(sat_az)
                        #alt.append(sat_alt)
                        ##t_rise.append(t_start)
                        #time_array.append(time_interp)
                        #sat_id.append(sat_num)
                        ##sat_id.extend(sat_ephem['sat_id'])
                        test = zip(sat_num, time_interp, sat_alt, sat_az)
                        #print(list(test))
                        
                        for e in list(test):
                            sat_num = e[0]
                            time_u = e[1]
                            altitiude = e[2]
                            azimuth = e[3]

                            if time_u >= unix_start and time_u <= unix_stop:
                                sat_name.append(sat_num)
                                time_unix.append(time_u)
                                alt.append(altitiude)
                                az.append(azimuth)



    
                    else:
                        pass
    #                break
    #        break                
    
    #_, time_array, alt, az, sat_id = zip(*sorted(zip(t_rise, time_array, alt, az, sat_id)))
    
    #chrono_ephem = {}
    #
    #chrono_ephem['sat_id']      = list(sat_id)
    #chrono_ephem['time_array']  = list(time_array)
    #chrono_ephem['sat_alt']     = list(alt)
    #chrono_ephem['sat_az']      = list(az)

    chrono_ephem = {}
    
    chrono_ephem['sat_id']      = list(sat_name)
    chrono_ephem['time_array']  = list(time_unix)
    chrono_ephem['sat_alt']     = list(alt)
    chrono_ephem['sat_az']      = list(az)


# Save collated sat ephem data to one big json file
with open(f'{json_dir}/{date}.json', 'w') as outfile:
    json.dump(chrono_ephem, outfile, indent=4)

## Save collated sat ephem data to one big json file
#with open(f'{json_dir}/ultimate_ephem_list.json', 'w') as outfile:
#    json.dump(chrono_ephem, outfile, indent=4)


