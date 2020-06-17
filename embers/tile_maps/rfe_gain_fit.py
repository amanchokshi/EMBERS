import sys
import json
import argparse
import matplotlib
import numpy as np

matplotlib.use("Agg")
from pathlib import Path
import scipy.optimize as opt
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from scipy.stats import median_absolute_deviation as mad
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.polynomial.polynomial as poly
from scipy.stats import binned_statistic
import scipy.optimize as opt

sys.path.append("../decode_rf_data")
from colormap import spectral

cmap = spectral()

parser = argparse.ArgumentParser(
    description="""
    Liner Fit for RF Explorer Gain variations
    """
)

parser.add_argument(
    "--rfe_dir",
    metavar="\b",
    default="../../outputs/tile_maps/rfe_gain/",
    help="Dir where RFE residual json files live. Default:../../outputs/tile_maps/rfe_gain/",
)

parser.add_argument(
    "--start_gain",
    metavar="\b",
    default=-50,
    type=int,
    help="Power at which RFE gain variations begin. Default: -50dBm",
)
parser.add_argument(
    "--stop_gain",
    metavar="\b",
    default=-30,
    type=int,
    help="Power at which RFE gain variations saturate. Default: -30dBm",
)

args = parser.parse_args()

start_gain = args.start_gain
stop_gain = args.stop_gain
rfe_dir = Path(args.rfe_dir)

# find all rfe_gain json files
gain_files = [item for item in rfe_dir.glob("*.json")]
names = [i.name.split(".")[0] for i in gain_files]

# Combine data from all RF Explorers
pass_data = []
pass_resi = []

for n, f in enumerate(gain_files):
    with open(f, "r") as data:
        rfe = json.load(data)
        pass_data.extend(rfe["pass_data"])
        pass_resi.extend(rfe["pass_resi"])

        # pass_data = rfe['pass_data']
        # pass_resi = rfe['pass_resi']
        #
        # plt.style.use('seaborn')
        # fig = plt.figure()
        # plt.scatter(pass_data, pass_resi, marker='.', alpha=0.7, color='seagreen')
        # plt.xlabel('Observed power [dBm]')
        # plt.ylabel('Residuals power [dB]')
        # plt.xlim([-80,-20])
        # plt.ylim([-20,20])
        # plt.tight_layout()
        # plt.savefig(f'../../outputs/paper_plots/{names[n]}.png', bbox_inches='tight')
        # plt.close()


fig = plt.figure()

plt.hexbin(pass_data, pass_resi, gridsize=121, cmap=cmap, alpha=0.99, zorder=0)

pass_data = np.array(pass_data)
pass_resi = np.array(pass_resi)

# Clean up huge noisy outliers
filtr = np.where(np.logical_and(pass_data <= -25, pass_data >= -65))
pass_data = pass_data[filtr]
pass_resi = pass_resi[filtr]

# Median of binned data
bin_med, bin_edges, binnumber = binned_statistic(
    pass_data, pass_resi, statistic="median", bins=16
)
bin_width = bin_edges[1] - bin_edges[0]
bin_centers = bin_edges[1:] - bin_width / 2

# Now look at data between -50, -30, where gains vary
filtr = np.where(np.logical_and(pass_data >= start_gain, pass_data <= stop_gain))
pass_data = pass_data[filtr]
pass_resi = pass_resi[filtr]

# Linear fit to RFE data between -50, -30 dBm
poly = np.polyfit(pass_data, pass_resi, 2)
# Mathematical function of fit, which can be evaluated anywhere
f = np.poly1d(poly)
# save polyfit to file
np.save(f"{rfe_dir}/rfe_gain_fit.npy", poly)

# x_f = [f.roots[0], -45, -40, -35, -30, -25, -20]
x_f = np.linspace(max(f.roots), -25, num=10)
y_f = f(x_f)

plt.plot(
    x_f,
    y_f,
    color="w",
    lw=2.1,
    marker="s",
    markeredgecolor="k",
    markeredgewidth=1.6,
    markersize=4.9,
    alpha=1,
    label="Gain fit",
    zorder=1,
)
plt.scatter(
    bin_centers,
    bin_med,
    marker="X",
    s=36,
    facecolors="#ee4540",
    lw=0.9,
    edgecolors="w",
    alpha=1,
    label="Median residuals",
    zorder=2,
)

leg = plt.legend(loc="lower right", frameon=True, markerscale=0.9, handlelength=1.4)
leg.get_frame().set_facecolor("#cccccc")
for l in leg.legendHandles:
    l.set_alpha(0.77)


plt.xlabel("Observed power [dBm]")
plt.ylabel("Residuals power [dB]")
plt.xlim([-65, -25])
plt.ylim([-10, 15])
plt.tick_params(axis="both", length=0)
plt.grid(color="#cccccc", alpha=0.36, lw=1.2)
plt.box(None)
plt.tight_layout()
plt.savefig(f"{rfe_dir}/rfe_gain_fit.png", bbox_inches="tight")
