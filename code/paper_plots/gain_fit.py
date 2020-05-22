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


plt.style.use('seaborn')

#plt.scatter(pass_data, pass_resi, marker='.', alpha=0.7, color='seagreen')
plt.hexbin(pass_data, pass_resi, gridsize=77, bins='log',cmap='Blues')

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

plt.plot(x_f, y_f, color='#ffba5a', lw=3)
plt.scatter(bin_centers, bin_med, color='#fe6845', marker='s', s=84)


plt.xlabel('Observed power [dB]')
plt.ylabel('Residuals [dB]')
plt.xlim([-40,-15])
plt.ylim([-20,15])
plt.tight_layout()
plt.savefig(f'global_gain_fit.png')


