import sys
import math
import json
import pytz
import json
import numpy as np
from pathlib import Path
from astropy.time import Time
from datetime import datetime, timedelta



def time_tree(start_date, stop_date, time_zone):
    "Start/end of 30 minuts obs in local, unix, gps"

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

    return (obs_time, obs_unix, obs_unix_end, obs_gps, obs_gps_end)


def read_json(out_dir, f_name):
    """Read ultimate_pointing_times.json"""

    with open(f'{out_dir}/{f_name}') as table:
        data = json.load(table)
        pointings   = data['grid_pt'] 
        start_gps   = data['start_gps']
        stop_gps    = data['stop_gps'] 
        obs_length  = data['obs_length']

    return (pointings, start_gps, stop_gps, obs_length)


def time_filter(s_rise, s_set, times):
    """How much of the 30 minute
    window was at a particular 
    pointing"""
    
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
    
    return occu


def obs_pointing(
        obs_time, 
        obs_gps, 
        obs_gps_end, 
        pointings,
        start_gps,
        stop_gps,
        start_date, 
        stop_date,
        out_dir
        ):

    """Find pointing of 30 min obs
    Classification of more than 70%
    of an obs is at one pointing"""
    
    point_0 = []
    point_2 = []
    point_4 = []
    point_41 = []
    
    # loop over each 30 min observation
    for i in range(len(obs_time)):
    
        # start, end times in gps format
        obs_window = [obs_gps[i], obs_gps_end[i]]
    
        # loop over mwa pointing observations from 
        # ultimate_pointing_times.json
        for j in range(len(start_gps)):
            occu = time_filter(start_gps[j], stop_gps[j], obs_window)
    
            if (occu != None and occu >= 0.7):
                if pointings[j] == 0:
                    point_0.append(obs_time[i])
                elif pointings[j] == 2:
                    point_2.append(obs_time[i])
                elif pointings[j] == 4:
                    point_4.append(obs_time[i])
                elif pointings[j] == 41:
                    point_41.append(obs_time[i])
                else:
                    pass


    # Create dictionary to be saved to json
    obs_pointings = {}
    
    obs_pointings['start_date'] = start_date
    obs_pointings['stop_date']  = stop_date
    obs_pointings['point_0']    = point_0
    obs_pointings['point_2']    = point_2
    obs_pointings['point_4']    = point_4
    obs_pointings['point_41']   = point_41
    
    
    with open(f'{out_dir}/obs_pointings.json', 'w') as outfile:
        json.dump(obs_pointings, outfile, indent=4)



if __name__ == "__main__":
    
    import argparse
    
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


    obs_time, obs_unix, obs_unix_end, obs_gps, obs_gps_end = time_tree(
            start_date,
            stop_date,
            time_zone)

    pointings, start_gps, stop_gps, obs_length = read_json(out_dir, f_name)

    obs_pointing(
        obs_time, 
        obs_gps, 
        obs_gps_end, 
        pointings,
        start_gps,
        stop_gps,
        start_date, 
        stop_date,
        out_dir
        )





