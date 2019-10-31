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

# open a Two Line Element (TLE) file
with open('TLE/21576.txt', 'r') as f:
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

# Midpoints between epochs of successive TLEs
epochs = np.asarray(epochs)
midpoints = (epochs[1:] + epochs[:-1]) / 2


# Create time ranges in which TLEs are optimally accurate
# for i in range(len(midpoints)):
#    if i == 0:
#        print(epochs[0], midpoints[0])
#    elif i < range(len(midpoints))[-1]:
#        print(midpoints[i-1], midpoints[i])
#    else:
#        print(midpoints[i], epochs[-1])


# find time between epochs/ midpoints in seconds, using Astropy Time
t1 = Time(epochs[0], format='jd').gps
t2 = Time(epochs[2], format='jd').gps
#t2 = Time(midpoints[0], format='jd').gps
dt = round(t1 - t2)

# Create a time array, every 20 seconds, starting at the first epoch, approximately
seconds = range(dt)
t = ts.utc(2019, 10, 28, 13, 00, seconds[::20])

# Define Satellite to be the first one in sats list
# Find possition of sat at each timestep of time array
satellite = sats[0]
orbit = (satellite - MWA).at(t)
alt, az, distance = orbit.altaz()

# Check if sat is above the horizon, return boolean array
above_horizon = alt.degrees >= 0

# Indicies of rare times that sats are above the horizon
indicies, = above_horizon.nonzero()

# Boundary times at which the sat either rises or sets
boundaries, = np.diff(above_horizon).nonzero()
# print(boundaries)

# Reshape into pairs rise & set indicies
passes = boundaries.reshape(len(boundaries) // 2, 2)
#print(passes)


def plot_sat(pass_indices):
    i, j = pass_indices

    # Set up the polar plot.
    plt.style.use('seaborn')
    ax = plt.subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,30,60,90], angle=0)
    #ax.set_xticklabels(['N', '', 'E', '', 'S', '', 'W', ''])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)

    # Draw line and labels.
    θ = az.radians
    r = alt.degrees
    ax.plot(θ[i:j+1], r[i:j+1], '-', linewidth=6, alpha=0.6, c='#1fab89')
    
    plt.tight_layout()
    plt.show()


plot_sat(passes[1])
