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

gain_files = [item for item in Path('./gain_fit/').glob('*.json')]

pass_data = []
pass_resi = []

for f in gain_files:
    with open(f, 'r') as data:
        rfe = json.load(data)
        pass_data.extend(rfe['pass_data'])
        pass_resi.extend(rfe['pass_resi'])


nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
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

#plt.style.use('seaborn')

fig = plt.figure(figsize=(3.6,2.4))

#plt.scatter(pass_data, pass_resi, marker='.', alpha=0.7, color='seagreen')
plt.hexbin(pass_data, pass_resi, gridsize=36, bins='log',cmap='Blues', alpha=0.9)

pass_data = np.array(pass_data)
pass_resi = np.array(pass_resi)
filtr = np.where(pass_data <= -20)
pass_data = pass_data[filtr]
pass_resi = pass_resi[filtr]

bin_med, bin_edges, binnumber = binned_statistic(pass_data, pass_resi, statistic='median', bins=11)
bin_width = (bin_edges[1] - bin_edges[0])
bin_centers = bin_edges[1:] - bin_width/2


z = np.polyfit(bin_centers, bin_med, 3)
f = np.poly1d(z)

x_f = np.linspace(-40, -15)
y_f = f(x_f)

plt.plot(x_f, y_f, color='#ffd800', lw=2.1, alpha=1, label='Gain fit')
plt.scatter(bin_centers, bin_med, color='#fe6845', marker='s', s=28, alpha=1, label='Gain Binned')

leg = plt.legend(loc="lower right", frameon=True, markerscale=1, handlelength=1)
leg.get_frame().set_facecolor('w')
for l in leg.legendHandles:
    l.set_alpha(1)


plt.xlabel('Observed power [dB]')
plt.ylabel('Residuals power [dB]')
plt.yticks([-20, -10, 0, 10])
plt.xlim([-40,-15])
plt.ylim([-20,15])
plt.tick_params(axis='both', length = 0)
plt.grid(color='w', alpha=0.42, lw=1.2)
plt.box(False)
plt.tight_layout()
plt.savefig(f'../../outputs/paper_plots/rfe_gain_fit.pdf', bbox_inches='tight')


