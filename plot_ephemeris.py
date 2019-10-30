import numpy as np
import skyfield as sf
from skyfield.api import Topos, load
from astropy.time import Time
import matplotlib.pyplot as plt

# Skyfield Timescale 
ts = load.timescale(builtin=True)

# Position of MWA site in Lat/Lon/Elevation
MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)

# open a Two Line Element (TLE) file
with open('TLE/21576.txt', 'r') as f:
    tle_list = [line.strip() for line in f.read().split('\n') if line is not '']

# List of satellite ephemeris from each set of TLE in the opened file
# Epoch: Time at which TLE is most accurate (in JD - Julian Date)
sats = []
epochs = []

for i in range(0,len(tle_list), 2):
    sat = sf.sgp4lib.EarthSatellite(tle_list[i],tle_list[i+1])
    epoch = sat.model.jdsatepoch
    sats.append(sat)
    epochs.append(epoch)

# Midpoints between epochs of successive TLEs
epochs = np.asarray(epochs)
midpoints = (epochs[1:] + epochs[:-1]) / 2


# Create time ranges in which TLEs are optimally accurate
#for i in range(len(midpoints)):
#    if i == 0:
#        print(epochs[0], midpoints[0])
#    elif i < range(len(midpoints))[-1]:
#        print(midpoints[i-1], midpoints[i])
#    else:
#        print(midpoints[i], epochs[-1])

t1 = Time(epochs[0], format='jd').gps    
t2 = Time(midpoints[0], format='jd').gps    
dt = round(t1 - t2)

seconds = range(dt)
t = ts.utc(2019, 10, 28, 13, 00, seconds)


satellite = sats[0]
orbit = (satellite - MWA).at(t)
alt, az, distance = orbit.altaz()
#print(alt)
#print(az)

above_horizon = alt.degrees > 0
#print(above_horizon)

indicies, = above_horizon.nonzero()
#print(indicies)

boundaries, = np.diff(above_horizon).nonzero()
#print(boundaries)


passes = boundaries.reshape(len(boundaries) // 2, 2)
#print(passes)

def plot_sky(pass_indices):
    i, j = pass_indices
    print('Rises:', t[i])
    print('Sets:', t[j])
       
    # Set up the polar plot.
    ax = plt.subplot(111, projection='polar')
    ax.set_rlim([0, 90])
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
                
    # Draw line and labels.
    θ = az.radians
    r = 90 - alt.degrees
    ax.plot(θ[i:j], r[i:j], 'ro--')
    #for k in range(i, j):
        #text = t[k].strftime('%H:%M')
        #ax.text(θ[k], r[k], text, ha='right', va='bottom')
    plt.show()

plot_sky(passes[0])
