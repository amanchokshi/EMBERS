import argparse
import concurrent.futures
from itertools import repeat
from pathlib import Path

import matplotlib
import numpy as np
from embers.tile_maps.beam_utils import (chisq_fit_gain,
                                         healpix_cardinal_slices, map_slices,
                                         poly_fit, rotate_map)
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.use("Agg")


parser = argparse.ArgumentParser(
    description="""
    Tile Slices paper plot
    """
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots/tile_slices/",
    help="Output directory. Default=../embers_out/paper_plots/tile_slices",
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
    "--ref_res",
    metavar="\b",
    default="../embers_out/tile_maps/null_test/ref_res.npy",
    help="Reference residuals from null test. default=../embers_out/tile_maps/null_test/ref_res.npy",
)
parser.add_argument(
    "--nside", metavar="\b", type=int, default=32, help="Healpix Nside. Default = 32",
)

args = parser.parse_args()

out_dir = Path(args.out_dir)
map_dir = Path(args.map_dir)
fee_map = args.fee_map
ref_res = args.ref_res
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


def plt_slice(
    fig=None,
    sub=(None, None, None),
    zen_angle=None,
    map_slice=None,
    map_error=None,
    model_slice=None,
    delta_pow=None,
    pow_fit=None,
    ref_res=None,
    slice_label=None,
    model_label=None,
    xlabel=False,
    ylabel=True,
    xlim=[-82, 82],
    ylim=[-26, 12],
    title=None,
):

    """Plot a slice of measured beam map with errorbars fit the fee beam model. Subplot with residual power.

    :param fig: Figure number
    :param sub: Subplot position, tuple of matplotlib indices. Ex: (1, 1, 1)
    :param zen_angle: Array of zenith angles
    :param map_slice: Array of beam powers from slice of beam map
    :param map_error: Array of errors on map_slice
    :param model_slice: Slice of beam model at given zenith angles
    :param delta_pow: Residual power between measured beam slice and model beam
    :param pow_fit: Polynomial fit to residual power
    :param slice_label: Label of measured beam slice
    :param model_label: Label of beam model slice
    :param xlabel: If True, plot X label of plot
    :param ylabel: If True, plot Y label of plot
    :param xlim: X limits on the plot. Default: [-82, 82]
    :param ylim: Y limits on the plot. Default: [-26, 12]
    :param title: Plot title

    :returns:
        - ax - :func:`~matplotlib.pyplot.subplot` object

    """

    ax = fig.add_subplot(sub[0], sub[1], sub[2])

    ax.errorbar(
        zen_angle,
        map_slice,
        yerr=map_error,
        fmt=".",
        color="#326765",
        ecolor="#7da87b",
        elinewidth=1.4,
        capsize=1.4,
        capthick=1.6,
        alpha=0.9,
        ms=7,
        label=slice_label,
    )

    ax.plot(
        zen_angle,
        model_slice,
        color="#c70039",
        linewidth=1.4,
        alpha=0.9,
        label=model_label,
    )

    ax.text(0.02, 0.88, title, horizontalalignment="left", transform=ax.transAxes)
    ax.set_xticks([-75, -50, -25, 0, 25, 50, 75])
    ax.set_yticks([0, -10, -20, -30, -40])

    leg = ax.legend(loc="lower center", frameon=True, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticklabels([])

    divider = make_axes_locatable(ax)
    dax = divider.append_axes("bottom", size="30%", pad=0.06)

    dax.scatter(zen_angle, delta_pow, marker=".", s=30, color="#27296d")
    dax.plot(zen_angle, pow_fit, linewidth=1.4, alpha=0.9, color="#ff8264")
    dax.fill_between(zen_angle, ref_res, alpha=0.9, color="#008891", edgecolor="black", zorder=0)
    dax.set_xticks([-75, -50, -25, 0, 25, 50, 75])
    dax.set_ylim([-10, 10])
    dax.set_xlim(xlim)
    #  dax.set_ylim([-5, 5])
    if ylabel is True:
        ax.set_ylabel("Power [dB]")
        dax.set_ylabel(r"$\Delta$ref [dB]")
    else:
        dax.set_yticklabels([])
        ax.set_yticklabels([])
    if xlabel is False:
        dax.set_xticklabels([])
    else:
        dax.set_xlabel("Zenith Angle [deg]")

    return ax


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


def plt_grid(tile_pair, fee_map, ref_res, nside, map_dir, out_dir):
    # This is an Awesome plot

    tile_name = tile_pair[0]
    ref_name = tile_pair[1]

    f_xx = f"{map_dir}/{tile_name}XX_{ref_name}XX_tile_maps.npz"
    f_yy = f"{map_dir}/{tile_name}YY_{ref_name}YY_tile_maps.npz"

    maps_xx = beam_slices(f_xx, fee_map, nside)
    maps_yy = beam_slices(f_yy, fee_map, nside)

    ref_res = np.load(ref_res, allow_pickle=True).flat[0]
    #  za = ref_res["za"]
    rf1_XX_NS = ref_res["rf1_XX_NS"]
    rf1_XX_EW = ref_res["rf1_XX_EW"]
    rf1_YY_NS = ref_res["rf1_YY_NS"]
    rf1_YY_EW = ref_res["rf1_YY_EW"]

    plt.style.use("seaborn")

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

    # ax1  = fig1.add_axes([0.01, 0.75, 0.29, 0.22])

    plt.subplot(4, 3, 1)
    plt.title(fr"($i$) NS slice of {tile_name}XX [zenith]", loc="left")
    plt_slice(
        fig=fig1,
        sub=[4, 3, 1],
        zen_angle=maps_xx[0][0][0][2],
        map_slice=maps_xx[0][0][2],
        map_error=maps_xx[0][0][0][1],
        model_slice=maps_xx[0][0][1][0],
        delta_pow=maps_xx[0][0][3],
        pow_fit=maps_xx[0][0][4],
        ref_res=rf1_XX_NS,
        slice_label="Tile NS",
        model_label="FEE NS",
        ylabel=True,
        xlim=[-92, 92],
        ylim=[-50, 5],
    )
    ax1 = plt.gca()
    ax1.set_xticks([-75, -50, -25, 0, 25, 50, 75])

    plt.subplot(4, 3, 2)
    plt.title(fr"($ii$) NS slice of {tile_name}XX [2]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 2],
        zen_angle=maps_xx[1][0][0][2],
        map_slice=maps_xx[1][0][2],
        map_error=maps_xx[1][0][0][1],
        model_slice=maps_xx[1][0][1][0],
        delta_pow=maps_xx[1][0][3],
        pow_fit=maps_xx[1][0][4],
        ref_res=rf1_XX_NS,
        slice_label="Tile NS",
        model_label="FEE NS",
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 3)
    plt.title(fr"($iii$) NS slice of {tile_name}XX [4]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 3],
        zen_angle=maps_xx[2][0][0][2],
        map_slice=maps_xx[2][0][2],
        map_error=maps_xx[2][0][0][1],
        model_slice=maps_xx[2][0][1][0],
        delta_pow=maps_xx[2][0][3],
        pow_fit=maps_xx[2][0][4],
        ref_res=rf1_XX_NS,
        slice_label="Tile NS",
        model_label="FEE NS",
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 4)
    plt.title(fr"($iv$) EW slice of {tile_name}XX [zenith]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 4],
        zen_angle=maps_xx[0][1][0][2],
        map_slice=maps_xx[0][1][2],
        map_error=maps_xx[0][1][0][1],
        model_slice=maps_xx[0][1][1][0],
        delta_pow=maps_xx[0][1][3],
        pow_fit=maps_xx[0][1][4],
        ref_res=rf1_XX_EW,
        slice_label="Tile EW",
        model_label="FEE EW",
        ylabel=True,
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 5)
    plt.title(fr"($v$) EW slice of {tile_name}XX [2]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 5],
        zen_angle=maps_xx[1][1][0][2],
        map_slice=maps_xx[1][1][2],
        map_error=maps_xx[1][1][0][1],
        model_slice=maps_xx[1][1][1][0],
        delta_pow=maps_xx[1][1][3],
        pow_fit=maps_xx[1][1][4],
        ref_res=rf1_XX_EW,
        slice_label="Tile EW",
        model_label="FEE EW",
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 6)
    plt.title(fr"($vi$) EW slice of {tile_name}XX [4]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 6],
        zen_angle=maps_xx[2][1][0][2],
        map_slice=maps_xx[2][1][2],
        map_error=maps_xx[2][1][0][1],
        model_slice=maps_xx[2][1][1][0],
        delta_pow=maps_xx[2][1][3],
        pow_fit=maps_xx[2][1][4],
        ref_res=rf1_XX_EW,
        slice_label="Tile EW",
        model_label="FEE EW",
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 7)
    plt.title(fr"($vii$) NS slice of {tile_name}YY [zenith]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 7],
        zen_angle=maps_yy[0][0][0][2],
        map_slice=maps_yy[0][0][2],
        map_error=maps_yy[0][0][0][1],
        model_slice=maps_yy[0][0][1][0],
        delta_pow=maps_yy[0][0][3],
        pow_fit=maps_yy[0][0][4],
        ref_res=rf1_YY_NS,
        slice_label="Tile NS",
        model_label="FEE NS",
        ylabel=True,
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 8)
    plt.title(fr"($viii$) NS slice of {tile_name}YY [2]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 8],
        zen_angle=maps_yy[1][0][0][2],
        map_slice=maps_yy[1][0][2],
        map_error=maps_yy[1][0][0][1],
        model_slice=maps_yy[1][0][1][0],
        delta_pow=maps_yy[1][0][3],
        pow_fit=maps_yy[1][0][4],
        ref_res=rf1_YY_NS,
        slice_label="Tile NS",
        model_label="FEE NS",
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 9)
    plt.title(fr"($ix$) NS slice of {tile_name}YY [4]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 9],
        zen_angle=maps_yy[2][0][0][2],
        map_slice=maps_yy[2][0][2],
        map_error=maps_yy[2][0][0][1],
        model_slice=maps_yy[2][0][1][0],
        delta_pow=maps_yy[2][0][3],
        pow_fit=maps_yy[2][0][4],
        ref_res=rf1_YY_NS,
        slice_label="Tile NS",
        model_label="FEE NS",
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 10)
    plt.title(
        fr"($x$) EW slice of {tile_name}YY [Zenith]", loc="left",
    )
    plt_slice(
        fig=fig1,
        sub=[4, 3, 10],
        zen_angle=maps_yy[0][1][0][2],
        map_slice=maps_yy[0][1][2],
        map_error=maps_yy[0][1][0][1],
        model_slice=maps_yy[0][1][1][0],
        delta_pow=maps_yy[0][1][3],
        pow_fit=maps_yy[0][1][4],
        ref_res=rf1_YY_EW,
        slice_label="Tile EW",
        model_label="FEE EW",
        ylabel=True,
        xlabel=True,
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 11)
    plt.title(fr"($xi$) EW slice of {tile_name}YY [2]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 11],
        zen_angle=maps_yy[1][1][0][2],
        map_slice=maps_yy[1][1][2],
        map_error=maps_yy[1][1][0][1],
        model_slice=maps_yy[1][1][1][0],
        delta_pow=maps_yy[1][1][3],
        pow_fit=maps_yy[1][1][4],
        ref_res=rf1_YY_EW,
        slice_label="Tile EW",
        model_label="FEE EW",
        xlabel=True,
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.subplot(4, 3, 12)
    plt.title(fr"($xii$) EW slice of {tile_name}YY [4]", loc="left"),
    plt_slice(
        fig=fig1,
        sub=[4, 3, 12],
        zen_angle=maps_yy[2][1][0][2],
        map_slice=maps_yy[2][1][2],
        map_error=maps_yy[2][1][0][1],
        model_slice=maps_yy[2][1][1][0],
        delta_pow=maps_yy[2][1][3],
        pow_fit=maps_yy[2][1][4],
        ref_res=rf1_YY_EW,
        slice_label="Tile EW",
        model_label="FEE EW",
        xlabel=True,
        xlim=[-92, 92],
        ylim=[-50, 5],
    )

    plt.tight_layout()
    plt.savefig(
        f"{out_dir}/{tile_name}_{ref_name}_map_slices.pdf", bbox_inches="tight", dpi=420
    )
    plt.close()


if __name__ == "__main__":

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(
            plt_grid,
            tile_pairs,
            repeat(fee_map),
            repeat(ref_res),
            repeat(nside),
            repeat(map_dir),
            repeat(out_dir),
        )
