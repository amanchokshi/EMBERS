import json
import argparse
import matplotlib
matplotlib.use
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=""" 
        Plot Satellite passes from json ephemeris files.
        """)
parser.add_argument('--sat', metavar='\b', help='Norad cat ID of satellite. Ex:23545')

args = parser.parse_args()
sat_num = args.sat

def sat_plot(alt, az):
    # Set up the polar plot.
    plt.style.use('seaborn')
    ax = plt.subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,10,20,30,40,50,60,70,80,90], angle=22)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.grid(color='darkslategrey', alpha=0.3)
    ax.plot(az, alt, '-', linewidth=1, alpha=0.4, c='seagreen') #1fab89
    plt.tight_layout()
    plt.show()        


with open('./outputs/ephem_json/{}.json'.format(sat_num), 'r') as ephem:
    sat_ephem = json.load(ephem)
    altitude = []
    azimuth = []
    
    for alt in sat_ephem['sat_alt']:
        altitude.extend(alt)
    for az in sat_ephem['sat_az']:
        azimuth.extend(az)
    sat_plot(altitude, azimuth)

        
