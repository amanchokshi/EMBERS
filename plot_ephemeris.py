import numpy as np
import skyfield as sf
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from astroplan.plots import plot_sky
from skyfield.api import Topos, load

# Skyfield Timescale
ts = load.timescale(builtin=True)

# Position of MWA site in Lat/Lon/Elevation
MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)


def load_tle(tle_path):
    '''Loads a tle file & returns list of sat objects and corresponding epochs''' 
    
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

    return (sats, epochs)


def epoch_ranges(epochs):
    '''Creates a time array with intervals corresponding to 
    epochs of best accuracy from the TLE file
    '''
    
    # Midpoints between epochs of successive TLEs
    epoch_list = np.asarray(epochs)
    midpoints = list((epoch_list[1:] + epoch_list[:-1]) / 2)
    epoch_range = [epochs[0]] + midpoints + [epochs[-1]]

    return epoch_range


def epoch_time_array(index_epoch, t_step):
    '''Create a skyfield time object [time vector] at which
    the sat position will be determined.
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
    seconds = range(dt)
    t_arr = ts.utc(int(year), int(month), int(day), int(hour), int(minute), seconds[::t_step])

    return (t_arr, index_epoch)



def sat_pass(t_arr, index_epoch):
    '''Calculates alt/az of sat at every position instant of time array (t_arr).
    Finds all the times that the sat is above the horizon, and returns the index
    t_arr at which the sat rose and set.'''

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


def plot_sat(pass_indices):
    '''Plots a satellite pass on a polar alt/az map'''
    i, j = pass_indices

    # Set up the polar plot.
    plt.style.use('seaborn')
    ax = plt.subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,30,60,90], angle=22)
    #ax.set_xticklabels(['N', '', 'E', '', 'S', '', 'W', ''])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    # Draw line and labels.
    θ = az.radians
    r = alt.degrees
    ax.plot(θ[i:j+1], r[i:j+1], '-', linewidth=6, alpha=0.6, c='#1fab89')
    
    plt.tight_layout()
    #plt.show()



sats, epochs = load_tle('TLE/21576.txt')
epoch_range = epoch_ranges(epochs)
t_arr, index_epoch = epoch_time_array(2, 20)
passes, alt, az = sat_pass(t_arr, index_epoch)


for i in range(len(passes)):
    if len(passes) == 0:
        print('No Satellite Passes')
    else:
        plot_sat(passes[i])

plt.show()
