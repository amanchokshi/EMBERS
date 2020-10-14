"""
SAVGOL INTERP.
--------------
"""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
from embers.rf_tools.align_data import savgol_interp
from embers.rf_tools.colormaps import spectral
from matplotlib import pyplot as plt

matplotlib.use("Agg")
_spec, _ = spectral()

parser = argparse.ArgumentParser(
    description="""
        Savgol Interpolation paper plot
        """
)

parser.add_argument(
    "--rf_dir",
    metavar="\b",
    default="../../tiles_data",
    help="Directory with raw rf data. Default=.../../tiles_data",
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots",
    help="Output Directory. Default=./embers_out/paper_plots",
)

args = parser.parse_args()
rf_dir = Path(args.rf_dir)
out_dir = Path(args.out_dir)

out_dir.mkdir(parents=True, exist_ok=True)

try:
    ch = 8
    (
        ref_ali,
        tile_ali,
        time_array,
        ref_power,
        tile_power,
        ref_time,
        tile_time,
    ) = savgol_interp(
        f"{rf_dir}/rf0XX/2019-09-15/rf0XX_2019-09-15-11:00.txt",
        f"{rf_dir}/S06XX/2019-09-15/S06XX_2019-09-15-11:00.txt",
        savgol_window_1=11,
        savgol_window_2=15,
        polyorder=2,
        interp_type="cubic",
        interp_freq=1,
    )

    plt.style.use("seaborn")

    nice_fonts = {
        # Use LaTeX to write all text
        # "text.usetex": True,
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

    fig = plt.figure(figsize=(3.6, 2.4))

    colors = _spec([0.14, 0.28])

    tile_t = tile_time - tile_time[0]
    time_array = time_array - time_array[0]

    med = np.median(tile_power)
    tile_p = tile_power - med
    tile_p_aligned = tile_ali - med

    plt.plot(
        time_array,
        tile_p_aligned[::, ch],
        linewidth=1,
        color=colors[0],
        #  color="#2c5d63",
        alpha=0.9,
        label="SavGol",
    )
    plt.scatter(
        tile_t,
        tile_p[::, ch],
        color=colors[1],
        #  color="#7fa998",
        marker=".",
        s=3,
        alpha=0.2,
        label="AUT raw",
    )

    leg = plt.legend(loc="upper right", frameon=True, markerscale=4, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.ylabel("Power [dB]")
    plt.xlabel("Time [s]")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/savgol.pdf", bbox_inches="tight")
    print(f"SAVGOL INTERP saved to {out_dir}")

except Exception as e:
    print(e)
    print("Missing input rf files. Check path to rf_dir")
