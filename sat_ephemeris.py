import numpy as np
import skyfield as sf
from astropy.time import Time
from skyfield.api import Topos, load


# Skyfield Timescale
ts = load.timescale(builtin=True)

# Position of MWA site in Lat/Lon/Elevation
MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)


def load_tle(tle_path):
    '''Loads a tle.
    
    Opens TLE file and uses skyfield to load each pair of TLE lines
    into an EarthSatellite object. Determines epoch for each.

    Agrs:
        the_path: Path to tle file

    Returns:
        sats: List of skyfield earth satellite objects
        epochs: Time at which each set of TLE lines is most accurate
    ''' 
    
    # open a Two Line Element (TLE) file
    with open(tle_path, 'r') as f:
        tle_list = [line.strip()
                    for line in f.read().split('\n') if line is not '']
    
    # List of satellite ephemeris from each set of TLE in the opened file
    # Epoch: Time at which TLE is most accurate (in JD - Julian Date)
    sats = []
    epochs = []
    
    for i in range(0, len(tle_list), 2):
        sat = sf.sgp4lib.EarthSatellite(tle_list[i], tle_list[i+1])
        epoch = sat.model.jdsatepoch
        sats.append(sat)
        epochs.append(epoch)

    sats = np.asarray(sats)
    epochs = np.asarray(epochs)

    return (sats, epochs)


def epoch_ranges(epochs):
    '''Creates a time array
    
    Time array with intervals corresponding to epochs of 
    best accuracy from the TLE pair of lines in TLE file.
     
     Args:
        epochs: List of epochs derived from TLE file

     Returns:
        epoch_range: List of times, between which TLE is most accurate
    '''
    
    # Midpoints between epochs of successive TLEs
    epoch_list = np.asarray(epochs)
    midpoints = list((epoch_list[1:] + epoch_list[:-1]) / 2)
    epoch_range = [epochs[0]] + midpoints + [epochs[-1]]
    epoch_list = np.asarray(epoch_range)
    
    return epoch_range


def epoch_time_array(index_epoch, cadence):
    '''Create a time array.
    
    Skyfield time object [time vector] at which
    the sat position will be determined.

    Args:
        index_epoch: Index of epoch range to be converted to time array
        t_step: Cadence of time array
        epoch_range: List of times, between which TLE is most accurate

    Returns:
        t_arr: Skyfield time object
        index_epoch: Same as above
    '''

    # find time between epochs/ midpoints in seconds, using Astropy Time
    t1 = Time(epoch_range[index_epoch], format='jd').gps
    t2 = Time(epoch_range[index_epoch+1], format='jd').gps
    dt = round(t1 - t2)
    
    t3 = Time(epoch_range[index_epoch], format='jd').iso
    date, time = t3.split()
    year, month, day = date.split('-')
    hour, minute,_ = time.split(':')
    
    # Create a time array, every 20 seconds, starting at the first epoch, approximately
    seconds = np.arange(0,dt,cadence)
    t_arr = ts.utc(int(year), int(month), int(day), int(hour), int(minute), seconds)

    return (t_arr, index_epoch)


def sat_pass(sats, t_arr, index_epoch):
    '''Finds satellite passes.
    
    Calculates alt/az of sat at every position instant of time 
    array (t_arr). Finds all the times that the sat is above
    the horizon, and returns the index t_arr at which the sat
    rose and set.
    
    Args:
        sats: list of Skyfield EarthSatellite objects
        t_arr: Skyfield time object
        index_epoch: Position in epoch_range
    
    Returns:
        passes: 2D array with pairs of indicies corresponding to rise/set of sat
        alt: Array of altitudes of sat at t_arr times
        az: Array of azimuths of sat at t_arr times
    '''

    # Define Satellite to be the first one in sats list
    # Find possition of sat at each timestep of time array
    satellite = sats[index_epoch]
    orbit = (satellite - MWA).at(t_arr)
    alt, az, _ = orbit.altaz()
    
    # Check if sat is above the horizon, return boolean array
    above_horizon = alt.degrees >= 0
    
    # Indicies of rare times that sats are above the horizon
    indicies, = above_horizon.nonzero()
    
    # Boundary times at which the sat either rises or sets
    boundaries, = np.diff(above_horizon).nonzero()
    
    if above_horizon[0] == True:
        boundaries = [indicies[0]]+ list(boundaries)
        boundaries = np.asarray(boundaries)
    
    if above_horizon[-1] == True:
        boundaries = list(boundaries) + [indicies[-1]]
        boundaries = np.asarray(boundaries)
    
    # Reshape into pairs rise & set indicies
    passes = boundaries.reshape(len(boundaries) // 2, 2)

    return (passes, alt, az)


def ephem_data(t_arr, pass_index, alt, az):
    '''Satellite Ephemeris Data.
    
    creates rise time, set times and alt, az arrays.
    
    Args:
        t_arr: Skyfield time array object
        pass_index: One pair of sat indicies from passed 2D array
        alt: Array of altitudes of sat at t_arr times
        az: Array of azimuths of sat at t_arr times

    Returns:
        t_rise: Rise time of sat in gps seconds (Astropy)
        t_set: Set time of sat in gps seconds (Astropy)
        sat_alt: Altitude of satellite while it is above the horizon. Array.
        sat_az: Azimuth of satellite while it is above the horizon. Array.
    '''

    i, j = pass_index
    t_rise = Time(t_arr[i].tt, format='jd').gps
    t_set = Time(t_arr[j].tt, format='jd').gps

    θ = list(az.radians)
    r = list(alt.degrees)

    sat_az = θ[i:j+1]
    sat_alt = r[i:j+1]
    
    return (t_rise, t_set, sat_alt, sat_az)


def sat_plot(sat_num, alt, az, num_passes):
    '''Plots satellite passes
    
    Args:
        alt: list of altitude values
        az: list of azimuth values
        num_passes: Number of satellite passes
    '''
    
    import matplotlib.pyplot as plt
    
    # Set up the polar plot.
    plt.style.use('seaborn')
    figure = plt.figure(figsize=(6,6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,10,20,30,40,50,60,70,80,90], angle=22)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_title('Satellite {} Sky Coverage: {} Passes'.format(sat_num,num_passes), y=1.08)
    ax.grid(color='#071011', linewidth=1.4, alpha=0.4)
    plt.tight_layout()
    
    for i in range(len(alt)):
        plt.plot(az[i], alt[i], '-', linewidth=1.4, alpha=0.4, color='#458323') #1fab89
    
    # Return plot for saving, showing, etc
    return plt


if __name__ == '__main__':
    
    import os
    import json
    import argparse

    parser = argparse.ArgumentParser(description="""
            Code which converts the TLE files downloaded with download_TLE.py
            into satellite ephemeris data: rise time, set time, alt/az arrays
            at a given time cadence. This is saved to a json file which will 
            be used to plot the satellite passes.
            """)
            
    parser.add_argument('--sat', metavar='\b', help='The Norad cat ID of satellite. Example: 21576')
    parser.add_argument('--tle', metavar='\b', default='./outputs/TLE', help='Directory where TLE files are saved. Default=./outputs/TLE')
    parser.add_argument('--cadence', metavar='\b', default=20, help='Rate at which sat alt/az is computed. Expensive! Default=20s')
    
    args = parser.parse_args()
    sat_name = args.sat
    tle_dir = args.tle
    cadence = int(args.cadence)
    
    tle_path = '{}/{}.txt'.format(tle_dir, sat_name)
    
    
    os.makedirs(os.path.dirname('./outputs/ephem_json/'), exist_ok=True)
    
    sat_ephem = {}
    sat_ephem['t_rise'] = []
    sat_ephem['t_set'] = []
    sat_ephem['sat_alt'] = []
    sat_ephem['sat_az'] = []
    
    sats, epochs = load_tle(tle_path)
 #   print(sats)
    epoch_range = epoch_ranges(epochs)
    for i in range(len(epoch_range) - 1):   
        t_arr, index_epoch = epoch_time_array(i, cadence)
        passes, alt, az = sat_pass(sats, t_arr, index_epoch) 
#        print(passes)
        for pass_index in passes:
            t_rise, t_set, sat_alt, sat_az = ephem_data(t_arr, pass_index, alt, az)
    
            sat_ephem['t_rise'].append(t_rise)
            sat_ephem['t_set'].append(t_set)
            sat_ephem['sat_alt'].append(sat_alt)
            sat_ephem['sat_az'].append(sat_az)
    
    with open('./outputs/ephem_json/{}.json'.format(sat_name), 'w') as outfile:
        json.dump(sat_ephem, outfile)    


