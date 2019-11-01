import json
import matplotlib
matplotlib.use
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def sat_plot(alt, az):
    # Set up the polar plot.
    plt.style.use('seaborn')
    ax = plt.subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,30,60,90], angle=22)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.plot(az, alt, '-', linewidth=1, alpha=0.2, c='#1fab89')
    plt.tight_layout()
    plt.show()        


with open('./outputs/ephem_json/23545.json', 'r') as ephem:
    sat_ephem = json.load(ephem)
    altitude = []
    azimuth = []
    
    for alt in sat_ephem['sat_alt']:
        altitude.extend(alt)
    for az in sat_ephem['sat_az']:
        azimuth.extend(az)
    sat_plot(altitude, azimuth)

        
