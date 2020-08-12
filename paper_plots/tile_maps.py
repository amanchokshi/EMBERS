"""
TILE MAPS.
----------
"""

import argparse
import concurrent.futures
from itertools import repeat
from pathlib import Path

import matplotlib
import numpy as np
from embers.rf_tools.colormaps import jade
from embers.tile_maps.beam_utils import plot_healpix
from matplotlib import pyplot as plt

matplotlib.use("Agg")
jade, _ = jade()

parser = argparse.ArgumentParser(
    description="""
    Tile Maps paper plot
    """
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots/tile_maps/",
    help="Output directory. Default=../embers_out/paper_plots/tile_maps",
)
parser.add_argument(
    "--map_dir",
    metavar="\b",
    default="../embers_out/tile_maps/tile_maps/tile_maps_clean/",
    help="Tile map directory. Default=../embers_out/tile_maps/tile_maps/tile_maps_clean",
)
parser.add_argument(
    "--fee_map",
    metavar="\b",
    default="../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
    help="Healpix FEE map of mwa tile. default=../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
)
parser.add_argument(
    "--nside", metavar="\b", type=int, default=32, help="Healpix Nside. Default = 32",
)

args = parser.parse_args()

out_dir = Path(args.out_dir)
map_dir = Path(args.map_dir)
fee_map = args.fee_map
nside = args.nside

# make output dir if it doesn't exist
out_dir.mkdir(parents=True, exist_ok=True)

tiles = [
    "S06",
    "S07",
    "S08",
    "S09",
    "S10",
    "S12",
    "S29",
    "S30",
    "S31",
    "S32",
    "S33",
    "S34",
    "S35",
    "S36",
]

refs = ["rf0", "rf1"]

tile_pairs = []
for r in refs:
    for t in tiles:
        tile_pairs.append([t, r])


def beam_maps(f):
    """Returns pairs of [AUT,FEE] maps, for each pointing"""
    f_name, _ = Path(f).name.split(".")
    t_name, r_name, _, _ = f_name.split("_")

    # pointings = ['0','2','4','41']
    pointings = ["0", "2", "4"]

    maps = []

    # load data from map .npz file
    tile_map = np.load(f, allow_pickle=True)
    fee_m = np.load(fee_map, allow_pickle=True)

    for p in pointings:

        tile = tile_map[p]

        if "XX" in t_name:
            fee = fee_m[p][0]
        else:
            fee = fee_m[p][1]

        # Visualize the tile map and diff map
        # healpix meadian map
        tile_med = np.asarray([(np.nanmedian(j) if j != [] else np.nan) for j in tile])

        residuals = tile_med - fee
        residuals[np.where(fee < -30)] = np.nan
        residuals[np.where(tile_med == np.nan)] = np.nan

        maps.append([tile_med, residuals])

    return maps


def plt_grid(tile_pair, map_dir, out_dir):
    # This is an Awesome plot

    tile_name, ref_name = tile_pair

    f_xx = f"{map_dir}/{tile_name}XX_{ref_name}XX_tile_maps.npz"
    f_yy = f"{map_dir}/{tile_name}YY_{ref_name}YY_tile_maps.npz"

    maps_xx = beam_maps(f_xx)
    maps_yy = beam_maps(f_yy)

    nice_fonts = {
        # Use LaTeX to write all text
        # "text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "font.size": 8,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }
    plt.rcParams.update(nice_fonts)

    fig1 = plt.figure(figsize=(7.6, 9.0))

    ax1 = fig1.add_axes([0.01, 0.75, 0.29, 0.22])
    plot_healpix(
        data_map=maps_xx[0][0],
        fig=fig1,
        title=fr"($i$) {tile_name}XX [Zenith]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )

    ax2 = fig1.add_axes([0.30, 0.75, 0.29, 0.22])
    plot_healpix(
        data_map=maps_xx[1][0],
        fig=fig1,
        title=fr"($ii$) {tile_name}XX [2]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )

    ax3 = fig1.add_axes([0.60, 0.75, 0.29, 0.22])
    plot_healpix(
        data_map=maps_xx[2][0],
        fig=fig1,
        title=fr"($iii$) {tile_name}XX [4]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )
    ax3 = plt.gca()
    image = ax3.get_images()[0]
    cax1 = fig1.add_axes([0.91, 0.75, 0.015, 0.22])
    fig1.colorbar(image, cax=cax1, label="Power [dB]")

    ax4 = fig1.add_axes([0.01, 0.50, 0.29, 0.22])
    plot_healpix(
        data_map=maps_xx[0][1],
        fig=fig1,
        title=fr"($iv$) {tile_name}XX - FEE [Zenith]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )

    ax5 = fig1.add_axes([0.30, 0.50, 0.29, 0.22])
    plot_healpix(
        data_map=maps_xx[1][1],
        fig=fig1,
        title=fr"($v$) {tile_name}XX - FEE [2]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )

    ax6 = fig1.add_axes([0.60, 0.50, 0.29, 0.22])
    plot_healpix(
        data_map=maps_xx[2][1],
        fig=fig1,
        title=fr"($vi$) {tile_name}XX - FEE [4]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )
    ax6 = plt.gca()
    image = ax6.get_images()[0]
    cax2 = fig1.add_axes([0.91, 0.50, 0.015, 0.22])
    fig1.colorbar(image, cax=cax2, label="Power [dB]")

    ax7 = fig1.add_axes([0.01, 0.25, 0.29, 0.22])
    plot_healpix(
        data_map=maps_yy[0][0],
        fig=fig1,
        title=fr"($vii$) {tile_name}YY [Zenith]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )

    ax8 = fig1.add_axes([0.30, 0.25, 0.29, 0.22])
    plot_healpix(
        data_map=maps_yy[1][0],
        fig=fig1,
        title=fr"($viii$) {tile_name}YY [2]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )

    ax9 = fig1.add_axes([0.60, 0.25, 0.29, 0.22])
    plot_healpix(
        data_map=maps_yy[2][0],
        fig=fig1,
        title=fr"($ix$) {tile_name}YY [4]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )
    ax9 = plt.gca()
    image = ax9.get_images()[0]
    cax3 = fig1.add_axes([0.91, 0.25, 0.015, 0.22])
    fig1.colorbar(image, cax=cax3, label="Power [dB]")

    ax10 = fig1.add_axes([0.01, 0.0, 0.29, 0.22])
    plot_healpix(
        data_map=maps_yy[0][1],
        fig=fig1,
        title=fr"($x$) {tile_name}YY - FEE [Zenith]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )

    ax11 = fig1.add_axes([0.30, 0.0, 0.29, 0.22])
    plot_healpix(
        data_map=maps_yy[1][1],
        fig=fig1,
        title=fr"($xi$) {tile_name}YY - FEE [2]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )

    ax12 = fig1.add_axes([0.60, 0.0, 0.29, 0.22])
    plot_healpix(
        data_map=maps_yy[2][1],
        fig=fig1,
        title=fr"($xii$) {tile_name}YY - FEE [4]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )
    ax12 = plt.gca()
    image = ax12.get_images()[0]
    cax4 = fig1.add_axes([0.91, 0.0, 0.015, 0.22])
    fig1.colorbar(image, cax=cax4, label="Power [dB]")

    plt.savefig(
        f"{out_dir}/{tile_name}_{ref_name}_maps.pdf", bbox_inches="tight", dpi=420
    )
    plt.close()


if __name__ == "__main__":

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(plt_grid, tile_pairs, repeat(map_dir), repeat(out_dir))
