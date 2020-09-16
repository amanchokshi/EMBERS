"""
WATERFALLS.
-----------

"""

import argparse
import time
from pathlib import Path
import pkg_resources

import matplotlib
import numpy as np
from embers.rf_tools import colormaps, rf_data
from matplotlib import pyplot as plt

matplotlib.use("Agg")

parser = argparse.ArgumentParser(
    description="""
        Waterfall paper plots
        """
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots",
    help="Output dir. Default=../embers_out/paper_plots",
)

args = parser.parse_args()
out_dir = Path(args.out_dir)

out_dir.mkdir(parents=True, exist_ok=True)

rf_path = pkg_resources.resource_filename("embers.kindle", "data/rf_data/rf0XX/2019-10-10/rf0XX_2019-10-10-02:30.txt")
tl_path = pkg_resources.resource_filename("embers.kindle", "data/rf_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt")

# read in raw data
ref_p, ref_t = rf_data.read_data(rf_file=rf_path)
tile_p, tile_t = rf_data.read_data(rf_file=tl_path)

# scale median to zero
ref_p_median = np.median(ref_p)
tile_p_median = np.median(tile_p)
r_image = ref_p - ref_p_median
t_image = tile_p - tile_p_median

# setting dynamic range of waterfall to be 30 dB above the median
vmin = 0
vmax = 30

# Custom spectral colormap
spec, _ = colormaps.spectral()

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
}

plt.rcParams.update(nice_fonts)

fig = plt.figure(figsize=(7, 4.2))

ax1 = fig.add_axes([0.00, 0.0, 0.46, 1])
ax2 = fig.add_axes([0.50, 0.0, 0.46, 1])
cax = fig.add_axes([0.98, 0.0, 0.015, 1])

ax1.imshow(t_image, vmin=vmin, vmax=vmax, interpolation="none", cmap=spec)
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
ax1.set_xlabel("Frequency [MHz]")

# Y-axis stuff
y_ax = t_image.shape[0]
y_ticks = np.arange(0, y_ax, t_step)
ax1.set_yticks(y_ticks)
ax1.set_yticklabels(t_tz)
ax1.set_ylabel("MWA local time [HH:MM]")

ax2.imshow(r_image, vmin=vmin, vmax=vmax, interpolation="none", cmap=spec)
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
ax2.set_xlabel("Frequency [MHz]")

# Y-axis stuff
ax2.set_yticklabels([])
ax2.set_yticks([])

plt.savefig(f"{out_dir}/waterfall.png", dpi=600, bbox_inches="tight")

print(f"WATERRFALL PLOT saved to {out_dir}")
