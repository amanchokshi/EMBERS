import argparse
import concurrent.futures
from itertools import repeat
from pathlib import Path

import matplotlib
import numpy as np
from embers.rf_tools.colormaps import spectral
from embers.tile_maps.beam_utils import (chisq_fit_gain,
                                         healpix_cardinal_slices, map_slices,
                                         poly_fit, rotate_map)
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.use("Agg")
_spec, _ = spectral()

parser = argparse.ArgumentParser(
    description="""
    Tile Slices paper plot
    """
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots/tile_slices/",
    help="Output directory. Default=../embers_out/paper_plots",
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

refs = ["rf1"]

tile_pairs = []
for r in refs:
    for t in tiles:
        tile_pairs.append([t, r])


def beam_slices(map_file, fee_map, nside):
    """Returns pairs of [NS,EW] slices of maps, for each pointing"""

    t_name, r_name, _, _ = Path(map_file).stem.split("_")

    pointings = ["0", "2", "4"]

    maps = []

    # load data from map .npz file
    tile_map = np.load(map_file, allow_pickle=True)
    fee_m = np.load(fee_map, allow_pickle=True)

    for p in pointings:

        tile = tile_map[p]

        if "XX" in t_name:
            fee = fee_m[p][0]
        else:
            fee = fee_m[p][1]

        # rotate maps so slices can be taken
        fee_r = rotate_map(nside, angle=-np.pi / 4, healpix_array=fee)
        tile_r = rotate_map(nside, angle=-np.pi / 4, healpix_array=tile)

        # slice the tile and fee maps along NS, EW
        # zenith angle thresh of 70 to determine fit gain factor
        NS_f, EW_f = healpix_cardinal_slices(nside, fee_r, 70)
        NS_t, EW_t = map_slices(nside, tile_r, 70)

        gain_NS = chisq_fit_gain(data=NS_t[0], model=NS_f[0])
        gain_EW = chisq_fit_gain(data=EW_t[0], model=EW_f[0])

        # slice the tile and fee maps along NS, EW.
        # the above gain factor is applied to full beam slices
        NS_fee, EW_fee = healpix_cardinal_slices(nside, fee_r, 90)
        NS_tile, EW_tile = map_slices(nside, tile_r, 90)

        # Scale the data so that it best fits the beam slice
        NS_tile_med = NS_tile[0] - gain_NS[0]
        EW_tile_med = EW_tile[0] - gain_EW[0]

        # delta powers
        del_NS = NS_tile_med - NS_fee[0]
        del_EW = EW_tile_med - EW_fee[0]

        # 3rd order poly fits for residuals
        fit_NS = poly_fit(NS_tile[2], del_NS, NS_tile[0], 3)
        fit_EW = poly_fit(EW_tile[2], del_EW, EW_tile[0], 3)

        maps.append(
            [
                [NS_tile, NS_fee, NS_tile_med, del_NS, fit_NS],
                [EW_tile, EW_fee, EW_tile_med, del_EW, fit_EW],
            ]
        )

    return maps


# Extract a bunch of data

NS_0 = {}
EW_0 = {}
NS_2 = {}
EW_2 = {}
NS_4 = {}
EW_4 = {}


for tile in tile_pairs:

    tile_name = tile[0]
    ref_name = tile[1]

    for pol in ["XX", "YY"]:
        f_tile = f"{map_dir}/{tile_name}{pol}_{ref_name}{pol}_tile_maps.npz"
        try:
            slices = beam_slices(f_tile, fee_map, nside)
            NS_0[f"{tile_name}{pol}"] = slices[0][0]
            EW_0[f"{tile_name}{pol}"] = slices[0][1]
            NS_2[f"{tile_name}{pol}"] = slices[1][0]
            EW_2[f"{tile_name}{pol}"] = slices[1][1]
            NS_4[f"{tile_name}{pol}"] = slices[2][0]
            EW_4[f"{tile_name}{pol}"] = slices[2][1]
        except Exception as e:
            print(e)


def slice_scatter(
    fig=None,
    sub=(None, None, None),
    data=None,
    pol=None,
    xlabel=False,
    ylabel=True,
    xlim=[-94, 94],
    ylim=[-55, 5],
    title=None,
):

    """Plot the scatter of various beam slices.

    :returns:
        - ax - :func:`~matplotlib.pyplot.subplot` object

    """

    ax = fig.add_subplot(sub[0], sub[1], sub[2])

    tile_keys = list(data.keys())
    colors = _spec(np.linspace(0.17, 0.9, len(tile_keys)))

    plt.plot(
        data[f"S06{pol}"][0][2],
        data[f"S06{pol}"][1][0],
        color="black",
        lw=1.4,
        alpha=0.9,
    )

    for i, k in enumerate(tile_keys):
        if pol in k:
            plt.scatter(
                data[k][0][2],
                data[k][2],
                label=k,
                s=14,
                color=colors[i],
                alpha=0.88,
                edgecolor="black",
                linewidth=0.2,
            )

    ax.text(0.02, 0.88, title, horizontalalignment="left", transform=ax.transAxes)
    ax.set_xticks([-75, -50, -25, 0, 25, 50, 75])
    ax.set_yticks([0, -10, -20, -30, -40, -50])

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if ylabel is True:
        ax.set_ylabel("Power [dB]")
    else:
        ax.set_yticklabels([])
    if xlabel is False:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel("Zenith Angle [deg]")

    return ax


# Plotting stuff
plt.style.use("seaborn")

nice_fonts = {
    "font.family": "sans-serif",
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "font.size": 8,
    "legend.fontsize": 6,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
}


plt.rcParams.update(nice_fonts)

fig1 = plt.figure(figsize=(7.6, 9.0))

plt.subplot(4, 3, 1)
plt.title(r"($i$) NS slice of XX beam [zenith]", loc="left")
slice_scatter(
    fig=fig1, sub=[4, 3, 1], data=NS_0, pol="XX",
)

plt.subplot(4, 3, 2)
plt.title(r"($ii$) NS slice of XX beam [2]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 2], data=NS_2, pol="XX", ylabel=False)

plt.subplot(4, 3, 3)
plt.title(r"($iii$) NS slice of XX beam [4]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 3], data=NS_4, pol="XX", ylabel=False)

plt.subplot(4, 3, 4)
plt.title(r"($iv$) EW slice of XX beam [zenith]", loc="left")
slice_scatter(
    fig=fig1, sub=[4, 3, 4], data=EW_0, pol="XX",
)

plt.subplot(4, 3, 5)
plt.title(r"($v$) EW slice of XX beam [2]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 5], data=EW_2, pol="XX", ylabel=False)


plt.subplot(4, 3, 6)
plt.title(r"($vi$) EW slice of XX beam [4]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 6], data=EW_4, pol="XX", ylabel=False)


plt.subplot(4, 3, 7)
plt.title(r"($vii$) NS slice of YY beam [zenith]", loc="left")
slice_scatter(
    fig=fig1, sub=[4, 3, 7], data=NS_0, pol="YY",
)

plt.subplot(4, 3, 8)
plt.title(r"($viii$) NS slice of YY beam [2]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 8], data=NS_2, pol="YY", ylabel=False)

plt.subplot(4, 3, 9)
plt.title(r"($ix$) NS slice of YY beam [4]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 9], data=NS_4, pol="YY", ylabel=False)

plt.subplot(4, 3, 10)
plt.title(r"($x$) EW slice of YY beam [zenith]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 10], data=EW_0, pol="YY", xlabel=True)

plt.subplot(4, 3, 11)
plt.title(r"($xi$) EW slice of YY beam [2]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 11], data=EW_2, pol="YY", ylabel=False, xlabel=True)


plt.subplot(4, 3, 12)
plt.title(r"($xii$) EW slice of YY beam [4]", loc="left")
slice_scatter(fig=fig1, sub=[4, 3, 12], data=EW_4, pol="YY", ylabel=False, xlabel=True)

plt.tight_layout()
plt.savefig(f"{out_dir}/beam_varience.pdf", bbox_inches="tight", dpi=420)
plt.close()
