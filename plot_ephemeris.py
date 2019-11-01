import json
import matplotlib
matplotlib.use
import matplotlib.pyplot as plt
plt.style.use('ggplot')


with open('./outputs/ephem_json/21576.json', 'r') as ephem:
    sat_ephem = json.load(ephem)
    #for i in sat_ephem['t_rise']:
    #    print(i)
    alt = sat_ephem['sat_alt'][0]
    az = sat_ephem['sat_az'][0]

# Set up the polar plot.
    plt.style.use('seaborn')
    ax = plt.subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0,30,60,90], angle=22)
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.plot(az, alt, '-', linewidth=2, alpha=0.6, c='#1fab89')
    plt.tight_layout()
    plt.show()        

#def plot_sat(pass_indices):
#    '''Plots a satellite pass on a polar alt/az map'''
#    i, j = pass_indices
#
#    # Set up the polar plot.
#    plt.style.use('seaborn')
#    ax = plt.subplot(111, polar=True)
#    ax.set_ylim(90, 0)
#    ax.set_rgrids([0,30,60,90], angle=22)
#    #ax.set_xticklabels(['N', '', 'E', '', 'S', '', 'W', ''])
#    ax.set_theta_zero_location('N')
#    ax.set_theta_direction(-1)
#
#    # Draw line and labels.
#    θ = az.radians
#    r = alt.degrees
#    ax.plot(θ[i:j+1], r[i:j+1], '-', linewidth=2, alpha=0.6, c='#1fab89')
#
#    plt.tight_layout()
#    #plt.show()        
