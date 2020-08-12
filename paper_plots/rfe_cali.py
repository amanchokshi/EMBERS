"""
RFE CALI.
---------
"""

import argparse
import json
from pathlib import Path

import numpy as np
from embers.rf_tools.colormaps import spectral
from matplotlib import pyplot as plt
from scipy.stats import binned_statistic

spec, _ = spectral()

parser = argparse.ArgumentParser(
    description="""
    RFE gain fit paper plot
    """
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

parser.add_argument(
    "--rfe_cali",
    metavar="\b",
    default="../embers_out/tile_maps/rfe_calibration",
    help="Directory where RF Explorer calibration data lives. Default=../embers_out/tile_maps/rfe_calibration",
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots",
    help="Output directory. Default=../embers_out/paper_plots",
)

args = parser.parse_args()
start_gain = args.start_gain
stop_gain = args.stop_gain
rfe_cali = args.rfe_cali
out_dir = args.out_dir


def rfe_collate_cali(start_gain, stop_gain, rfe_cali_dir):
    """Collate RF Explorer gain calibration data from all MWA tile pairs, and plot a gain solution.

    :param start_gain: Power at which RFE gain variations begin. Ex: -50dBm
    :param stop_gain: Power at which RFE gain variations saturate. Ex: -30dBm
    :param rfe_cali_dir: Path to directory which contains gain calibration data saved by :func:`~embers.tile_maps.tile_maps.rfe_calibration`

    :returns:
        - Plot of global gain calibration solution and polynomial fit saved to :samp:`.npz` in the out_dir

    """

    # find all rfe_gain json files
    gain_files = [item for item in Path(rfe_cali_dir).glob("*.json")]

    # Combine data from all RF Explorers
    pass_data = []
    pass_resi = []

    for n, f in enumerate(gain_files):
        with open(f, "r") as data:
            rfe = json.load(data)
            pass_data.extend(rfe["pass_data"])
            pass_resi.extend(rfe["pass_resi"])

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

    plt.figure(figsize=(3.6, 2.4))

    plt.hexbin(pass_data, pass_resi, gridsize=121, cmap=spec, alpha=0.99, zorder=0)

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

    # x_f = [f.roots[0], -45, -40, -35, -30, -25, -20]
    x_f = np.linspace(max(f.roots), -25, num=10)
    y_f = f(x_f)

    plt.plot(
        x_f,
        y_f,
        color="w",
        lw=1.8,
        marker="s",
        markeredgecolor="k",
        markeredgewidth=1.4,
        markersize=4.2,
        alpha=1,
        label="Gain fit",
        zorder=1,
    )
    plt.scatter(
        bin_centers,
        bin_med,
        marker="X",
        s=24,
        facecolors="#ee4540",
        lw=0.7,
        edgecolors="w",
        alpha=1,
        label="Median residuals",
        zorder=2,
    )

    leg = plt.legend(loc="lower right", frameon=True, markerscale=0.96, handlelength=1.4)
    leg.get_frame().set_facecolor("#cccccc")
    for le in leg.legendHandles:
        le.set_alpha(0.9)

    plt.xlabel("Observed power [dBm]")
    plt.ylabel("Residuals power [dB]")
    plt.xlim([-65, -25])
    plt.ylim([-10, 15])
    plt.tick_params(axis="both", length=0)
    plt.grid(color="#cccccc", alpha=0.36, lw=1.2)
    plt.box(None)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/rfe_gain_fit.pdf", bbox_inches="tight")


rfe_collate_cali(start_gain, stop_gain, rfe_cali)
print(f"RFE CALI saved to {out_dir}")
