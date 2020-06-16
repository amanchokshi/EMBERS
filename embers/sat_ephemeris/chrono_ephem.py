import sys
import math
import json
import pytz
import numpy as np

from pathlib import Path
from scipy import interpolate
from datetime import datetime, timedelta


def obs_times(time_zone, start_date, stop_date):
    """Time stuff for 30 minute observations

    returns:
        obs_time        : YYYY-MM-DD-HH:MM of start of 30 min observations
        obs_unix        : Start of 30 min obs in unix time
        obs_unix_end    : End of 30 min obs in unix time
    """

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
    
    return (obs_time, obs_unix, obs_unix_end)


def write_json(data, filename=None, out_dir=None):
    '''writes data to json file in output dir'''

    with open(f'{out_dir}/{filename}','w') as f: 
        json.dump(data, f, indent=4) 


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
            
    # Create interpolation functions. Math functions, not Python!
    alt_interp = interpolate.interp1d(t_array, s_alt, kind=interp_type)
    
    # The next step was a bit tricky. Azimuth may wrap around.  
    # This may lead to a discontinuity between 0, 2π
    # We deal with this, by unwrapping the angles
    
    # This extends the angles beyond 2π, if the angles cross the discontinuity
    s_az_cont = np.unwrap(s_az)
    az_interp  = interpolate.interp1d(t_array, s_az_cont, kind=interp_type)
    
    # Makes start and end times clean integers
    # Also ensures that the interp range is inclusive of data points
    start = math.ceil(t_array[0])
    stop = math.floor(t_array[-1])
    
    # Create time array, at which to evaluate alt/az of sat
    time_interp = list(np.double(np.arange(start, stop, (1/interp_freq))))
    
    sat_alt = list(alt_interp(time_interp))
    
    # The modulus division by 2π, un-does the np.unwrap
    sat_az = list(az_interp(time_interp)%(2*np.pi))
    
    return(time_interp, sat_alt, sat_az)




if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="""
            Collates satellite pass data from all ephem json files.
            The native skyfiled gps timestamps are converted to unix
            timestamps to match the output of the rf explorers. The 
            alt, az data is interpolated to match the cadence of 
            align_data.py. Make a json file with all the passes 
            from each 30 min observation. This will help in the 
            next stage, where we identify all sats in each obs.
            """)
    
    parser.add_argument('--json_dir', metavar='\b', default='./../../outputs/sat_ephemeris/ephem_json/', help='Directory where ephem json files live. Default=./../../outputs/sat_ephemeris/ephem_json/')
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/sat_ephemeris/chrono_json/', help='Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/')
    parser.add_argument('--interp_type', metavar='\b', default='cubic', help='Type of interpolation. Ex: Cubic,Linear. Default=cubic')
    parser.add_argument('--interp_freq', metavar='\b', type=int, default=1, help='Frequency at which to interpolate, in Hertz. Must be the same as used in align_data.py. Default=1')
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
    
    
    # creates output dir, if it doesn't exist
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    
    # Save log file 
    sys.stdout = open(f'{out_dir}/logs_{start_date}_{stop_date}.txt', 'a')

    obs_time, obs_unix, obs_unix_end = obs_times(time_zone, start_date, stop_date)
    

    # Lets make the a json file for each 30 min observation, with an empty list
    data = []
    for i in range(len(obs_time)):
        write_json(data, filename=f'{obs_time[i]}.json', out_dir=out_dir)
    
    
    # Finds all sat ephem json files, and loops over them
    for json_path in list(Path(json_dir).glob('*.json')):
        
        with open(json_path) as ephem:
            print(json_path)
            sat_ephem = json.load(ephem)
            
            # Extract data from json dictionary
            t_array = sat_ephem['time_array']
            s_alt   = sat_ephem['sat_alt']
            s_az    = sat_ephem['sat_az']
            s_id  = sat_ephem['sat_id']
        
            # here, we're looping over each satellite pass with a single sat ephem file 
            # to check which observation window it falls in
            
            for pass_idx in range(len(t_array)):
                
                # Don't consider passes with less than 4 samples (~ 1 min = 3*20s)
                if len(t_array[pass_idx]) > 3:
                    time_interp, sat_alt, sat_az = interp_ephem(
                            t_array[pass_idx],
                            s_alt[pass_idx],
                            s_az[pass_idx],
                            interp_type,
                            interp_freq)
                    
                    # Loop over every sat pass, within the sat ephem file
                    for obs_int in range(len(obs_unix)):
                        
                        sat_ephem = {}
                        sat_ephem['sat_id'] = [s_id[0]]
                        sat_ephem['time_array'] = []
                        sat_ephem['sat_alt'] = []
                        sat_ephem['sat_az'] = []
                    
    
                        # Case I: Satpass occurs completely within the 30min observation
                        if (obs_unix[obs_int] < time_interp[0] and
                                obs_unix_end[obs_int] > time_interp[-1]):
                                
                            # append the whole pass to the dict
                            sat_ephem['time_array'].extend(time_interp)
                            sat_ephem['sat_alt'].extend(sat_alt)
                            sat_ephem['sat_az'].extend(sat_az)
                            #print(f'{pass_idx}: I.   {obs_time[obs_int]}')
                    
                    
                        # Case II: Satpass begins before the obs, but ends within it
                        elif (obs_unix[obs_int] > time_interp[0] and
                                obs_unix[obs_int] < time_interp[-1] and
                                obs_unix_end[obs_int] > time_interp[-1]):
                                
                            # find index of time_interp == obs_unix
                            start_idx = (np.where(np.asarray(time_interp) == obs_unix[obs_int]))[0][0]
                            
                            # append the end of the pass which is within the obs
                            sat_ephem['time_array'].extend(time_interp[start_idx:])
                            sat_ephem['sat_alt'].extend(sat_alt[start_idx:])
                            sat_ephem['sat_az'].extend(sat_az[start_idx:])
                    
                            #print(f'{pass_idx}: II.  {obs_time[obs_int]}')
                    
                        # Case III: Satpass begins within the obs, but ends after it
                        elif (obs_unix_end[obs_int] > time_interp[0] and 
                                obs_unix_end[obs_int] < time_interp[-1] and 
                                obs_unix[obs_int] < time_interp[0]):
                            
                            # find index of time_interp == obs_unix_end
                            stop_idx = (np.where(np.asarray(time_interp) == obs_unix_end[obs_int]))[0][0]
                            
                            # append the end of the pass which is within the obs
                            sat_ephem['time_array'].extend(time_interp[:stop_idx+1])
                            sat_ephem['sat_alt'].extend(sat_alt[:stop_idx+1])
                            sat_ephem['sat_az'].extend(sat_az[:stop_idx+1])
                            
                            #print(f'{pass_idx}: III. {obs_time[obs_int]}')
                    
                        # doesn't create json if there are no satellite passes within it
                        if sat_ephem['time_array'] != []:
                            
                            print(f'Satellite {s_id[0]} in {obs_time[obs_int]}')
                             
                            # open the relevant json file and loads contents to 'data_json'
                            with open(f'{out_dir}/{obs_time[obs_int]}.json') as json_file:
                                data_json = json.load(json_file)
                                
                                # append new satpass ephem data to data_json
                                data_json.append(sat_ephem)
    
                                # write the combined data back to the original file
                                write_json(data_json, filename=f'{obs_time[obs_int]}.json', out_dir=out_dir)
                                
                                # clear data_json
                                data_json = []


