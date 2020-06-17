import os
import sys
import time
import argparse
import numpy as np

sys.path.append("../../decode_rf_data")
import matplotlib.pyplot as plt
from colormap import spectral
import rf_data as rf

parser = argparse.ArgumentParser(
    description="""
        Will plot waterfall plots
        """
)

parser.add_argument(
    "--rf_dir",
    metavar="\b",
    default="./../../../data/",
    help="Path to rf data directory. Default=./../../../data/",
)
parser.add_argument(
    "--ref_name",
    metavar="\b",
    default="rf0XX_2019-10-10-02:30",
    help="Name of ref data file. Default=rf0XX_2019-10-10-02:30",
)
parser.add_argument(
    "--tile_name",
    metavar="\b",
    default="S10XX_2019-10-10-02:30",
    help="Name of tile data file. Default=S10XX_2019-10-10-02:30",
)
parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./../../outputs/paper_plots/",
    help="Output dir. Default=./../../outputs/paper_plots/",
)

args = parser.parse_args()
rf_dir = args.rf_dir
ref_name = args.ref_name
tile_name = args.tile_name
out_dir = args.out_dir

# Make output dir if it doesn't exist
os.makedirs(os.path.dirname(out_dir), exist_ok=True)


# read in raw data
ref_p, ref_t = rf.read_data(f"{rf_dir}/{ref_name}.txt")
tile_p, tile_t = rf.read_data(f"{rf_dir}/{tile_name}.txt")

# scale median to zero
ref_p_median = np.median(ref_p)
tile_p_median = np.median(tile_p)
r_image = ref_p - ref_p_median
t_image = tile_p - tile_p_median

# setting dynamic range of waterfall to be 30 dB above the median
vmin = 0
vmax = 30

# Custom spectral colormap
cmap = spectral()

nice_fonts = {
    # Use LaTeX to write all text
    "text.usetex": True,
    "font.family": "sans-serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "ytick.color": "#696969",
    "xtick.color": "#696969",
    "axes.labelcolor": "#696969",
    "axes.edgecolor": "#696969",
}

plt.rcParams.update(nice_fonts)

fig = plt.figure(figsize=(7, 4.2))
# fig, axs = plt.subplots(1,2, figsize=(7,5))
# axs = axs.ravel()

ax1 = fig.add_axes([0.00, 0.0, 0.46, 1])
ax2 = fig.add_axes([0.50, 0.0, 0.46, 1])
cax = fig.add_axes([0.98, 0.0, 0.015, 1])

ax1.imshow(t_image, vmin=vmin, vmax=vmax, interpolation="none", cmap=cmap)
ax1.set_aspect("auto")
image = ax1.get_images()[0]
cbar = fig.colorbar(image, cax=cax, label="Power [dB]")

# Number of time steps on y-axis
number_t = 5
t_step = int(len(tile_t) / (number_t - 1))
times = list(tile_t)
times = times[::t_step]

t_tz = []

# Convert UNIX time to local HH:MM time
for i in range(len(times)):

    perth_t = float(times[i]) + 28800  # 28800=+8GMT @ PERTH
    hms = time.strftime("%H:%M", time.gmtime(perth_t))
    t_tz.append(hms)

# Frequency: x-axis
start_freq = 137.15
stop_freq = 138.55

# X-axis stuff
x_ax = t_image.shape[1]
freqs = np.arange(start_freq, stop_freq, 0.25)
x_ticks = np.arange(0, x_ax, (0.25 / 0.0125))  # .0125MHz/ch
ax1.set_xticks(x_ticks)
ax1.set_xticklabels(freqs)
ax1.set_xlabel("Freqency [MHz]")

# Y-axis stuff
y_ax = t_image.shape[0]
y_ticks = np.arange(0, y_ax, t_step)
ax1.set_yticks(y_ticks)
ax1.set_yticklabels(t_tz)
ax1.set_ylabel("MWA local time [HH:MM]")


ax2.imshow(r_image, vmin=vmin, vmax=vmax, interpolation="none", cmap=cmap)
ax2.set_aspect("auto")
# Number of time steps on y-axis
number_t = 5
t_step = int(len(ref_t) / (number_t - 1))
times = list(ref_t)
times = times[::t_step]

t_tz = []

# Convert UNIX time to local HH:MM time
for i in range(len(times)):

    perth_t = float(times[i]) + 28800  # 28800=+8GMT @ PERTH
    hms = time.strftime("%H:%M", time.gmtime(perth_t))
    t_tz.append(hms)

# Frequency: x-axis
start_freq = 137.15
stop_freq = 138.55

# X-axis stuff
x_ax = r_image.shape[1]
freqs = np.arange(start_freq, stop_freq, 0.25)
x_ticks = np.arange(0, x_ax, (0.25 / 0.0125))  # .0125MHz/ch
ax2.set_xticks(x_ticks)
ax2.set_xticklabels(freqs)
ax2.set_xlabel("Freqency [MHz]")

# Y-axis stuff
ax2.set_yticklabels([])
ax2.set_yticks([])


plt.savefig(f"waterfall.png", dpi=144, transparent=True, bbox_inches="tight")
