import sys
import json
import argparse
import matplotlib
import numpy as np
matplotlib.use('Agg')
from pathlib import Path
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from scipy.stats import median_absolute_deviation as mad
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.polynomial.polynomial as poly
from scipy.stats import binned_statistic
import scipy.optimize as opt

sys.path.append('../decode_rf_data')
from colormap import spectral
cmap = spectral()

gain_files = [item for item in Path('../../outputs/paper_plots/gain_fit/').glob('*.json')]
names = [i.name.split('.')[0] for i in gain_files]

pass_data = []
pass_resi = []

for n,f in enumerate(gain_files):
    with open(f, 'r') as data:
        rfe = json.load(data)
        pass_data.extend(rfe['pass_data'])
        pass_resi.extend(rfe['pass_resi'])
        
nice_fonts = {
        # Use LaTeX to write all text
        #"text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 10,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        }

plt.rcParams.update(nice_fonts)

fig = plt.figure(figsize=(3.6,2.4))
#fig = plt.figure()

plt.hexbin(pass_data, pass_resi, gridsize=77,cmap=cmap, alpha=0.99, zorder=0)

pass_data = np.array(pass_data)
pass_resi = np.array(pass_resi)
filtr = np.where(np.logical_and(pass_data <= -20, pass_data >=-70))
pass_data = pass_data[filtr]
pass_resi = pass_resi[filtr]

# Median of binned data
bin_med, bin_edges, binnumber = binned_statistic(pass_data, pass_resi, statistic='median', bins=20)
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2

# Now look at data between -50, -30, where gains vary
filtr = np.where(np.logical_and(pass_data <= -28, pass_data >=-52))
pass_data = pass_data[filtr]
pass_resi = pass_resi[filtr]

z = np.polyfit(pass_data, pass_resi, 1)
f = np.poly1d(z)

x_f = [f.roots[0], -45, -40, -35, -30, -25, -20]
y_f = f(x_f)

#plt.plot(x_f, y_f, color='w', lw=2.1, marker='X', markeredgecolor='k', markeredgewidth=0.7, markersize=2.1, alpha=1, label='Gain fit', zorder=1)
plt.plot(x_f, y_f, color='w', lw=1.8, alpha=1, label='Gain fit', zorder=1)
plt.scatter(bin_centers, bin_med, marker='h', s=14,  facecolors='w',lw=0.7, edgecolors='k', alpha=1, label='Median residuals', zorder=2)
plt.scatter(f.roots[0], 0, marker='o', s=9,  facecolors='w', alpha=1)

leg = plt.legend(loc="lower right", frameon=True,framealpha=0.7, markerscale=0.9, handlelength=1)
leg.get_frame().set_facecolor('#cccccc')
for l in leg.legendHandles:
    l.set_alpha(1)


plt.xlabel('Observed power [dBm]')
plt.ylabel('Residuals power [dB]')
plt.xlim([-70,-20])
plt.ylim([-10,15])
plt.tick_params(axis='both', length = 0)
plt.grid(color='#cccccc', alpha=0.21, lw=1)
plt.box(None)
plt.tight_layout()
plt.savefig(f'../../outputs/paper_plots/rfe_gain_fit.pdf', bbox_inches='tight')

