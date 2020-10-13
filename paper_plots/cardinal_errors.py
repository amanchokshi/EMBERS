# Script to determine the differece in power between the FEE/measured maps along cardinal slices
# More power measured alone more sensitive axis of the dipole - perpendicular to dipole
# And vice-versa

# Rotation of reference tiles

import argparse
from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
from astropy.stats import mad_std as mad
from embers.rf_tools.colormaps import jade, spectral
from embers.tile_maps.beam_utils import (chisq_fit_gain,
                                         healpix_cardinal_slices, map_slices,
                                         poly_fit, rotate_map)
from healpy.rotator import euler_matrix_new as euler
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

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

        #  fee_r[PB_0] = np.nan
        #  tile_r[PB_0] = np.nan

        # slice the tile and fee maps along NS, EW
        # zenith angle thresh of 70 to determine fit gain factor
        NS_f, EW_f = healpix_cardinal_slices(nside, fee_r, 70)
        NS_t, EW_t = map_slices(nside, tile_r, 70)

        gain_NS = chisq_fit_gain(data=NS_t[0], model=NS_f[0])
        gain_EW = chisq_fit_gain(data=EW_t[0], model=EW_f[0])

        # slice the tile and fee maps along NS, EW.
        # the above gain factor is applied to full beam slices
        NS_fee, EW_fee = healpix_cardinal_slices(nside, fee_r, 85)
        NS_tile, EW_tile = map_slices(nside, tile_r, 85)

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


def reproject_map(nside, phi, healpix_array=None):
    """Reproject beam map to be centered around healpix pixel instead of zenith.

    :param nside: Healpix nside
    :param phi: Angle to rotate beam by, from E --> N (anticlockwise)
    :param healpix_array: Input healpix array to be rotated
    :returns:
        - Reprojected healpix map
    """

    vec = hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))
    eu_mat = euler(-phi, 0, 0, deg=True)
    rot_map = hp.rotator.rotateVector(eu_mat, vec)
    new_hp_inds = hp.vec2pix(nside, rot_map[0], rot_map[1], rot_map[2])

    return healpix_array[new_hp_inds]


def slice_residual():

    NS_r_XX = []
    NS_r_YY = []
    EW_r_XX = []
    EW_r_YY = []

    for tile in tiles:

        for pol in ["XX", "YY"]:

            f_tile = f"{map_dir}/{tile}{pol}_rf1{pol}_tile_maps.npz"

            try:
                map_slices = beam_slices(f_tile, fee_map, nside)

                #  pointings = ["0", "2", "4"]
                pointings = ["0"]

                for i, p in enumerate(pointings):
                    slices = map_slices[i]
                    NS_slices = slices[0]
                    fee_slice = NS_slices[1][0]
                    NS_med = NS_slices[0][0]
                    nulls = np.where(fee_slice < -30)
                    fee_slice[nulls] = np.nan
                    NS_med[nulls] = np.nan
                    NS_resi = np.array(NS_med - fee_slice)
                    if pol == "XX":
                        NS_r_XX.append(NS_resi)
                    else:
                        NS_r_YY.append(NS_resi)

                    EW_slices = slices[1]
                    fee_slice = EW_slices[1][0]
                    EW_med = EW_slices[0][0]
                    nulls = np.where(fee_slice < -30)
                    fee_slice[nulls] = np.nan
                    EW_med[nulls] = np.nan
                    EW_resi = np.array(EW_med - fee_slice)
                    if pol == "XX":
                        EW_r_XX.append(EW_resi)
                    else:
                        EW_r_YY.append(EW_resi)

            except Exception as e:
                print(e)

    NS_XX_res = np.nanmedian(np.array(NS_r_XX), axis=0)
    NS_YY_res = np.nanmedian(np.array(NS_r_YY), axis=0)
    EW_XX_res = np.nanmedian(np.array(EW_r_XX), axis=0)
    EW_YY_res = np.nanmedian(np.array(EW_r_YY), axis=0)

    NS_XX_mad = mad(np.array(NS_r_XX), axis=0)
    NS_YY_mad = mad(np.array(NS_r_YY), axis=0)
    EW_XX_mad = mad(np.array(EW_r_XX), axis=0)
    EW_YY_mad = mad(np.array(EW_r_YY), axis=0)

    # Only used to extract za arrays
    f_tile = f"{map_dir}/S06XX_rf0XX_tile_maps.npz"
    map_slices = beam_slices(f_tile, fee_map, nside)

    za = map_slices[0][0][0][2]

    NS_XX_fit = poly_fit(za, NS_XX_res, NS_XX_res, 2)
    NS_YY_fit = poly_fit(za, NS_YY_res, NS_YY_res, 2)
    EW_XX_fit = poly_fit(za, EW_XX_res, EW_XX_res, 2)
    EW_YY_fit = poly_fit(za, EW_YY_res, EW_YY_res, 2)

    # Plotting stuff
    #  plt.style.use("seaborn")

    #  nice_fonts = {
    #  "font.family": "sans-serif",
    #  "axes.labelsize": 8,
    #  "axes.titlesize": 9,
    #  "font.size": 8,
    #  "legend.fontsize": 6,
    #  "xtick.labelsize": 8,
    #  "ytick.labelsize": 8,
    #  }

    #  plt.rcParams.update(nice_fonts)

    #  plt.figure(figsize=(3.6, 2.4))

    #  colors = _spec(np.linspace(0.17, 0.9, 4))
    colors = _spec([0.14, 0.77, 0.66, 0.35])

    #  plt.figure(figsize=(6, 5))
    #  plt.errorbar(
    #  za,
    #  NS_XX_res,
    #  yerr=NS_XX_mad,
    #  fmt=".",
    #  ms=7,
    #  alpha=0.88,
    #  color=colors[0],
    #  elinewidth=1.4,
    #  capsize=1.4,
    #  capthick=1.6,
    #  )
    #  plt.errorbar(
    #  za,
    #  NS_YY_res,
    #  yerr=NS_YY_mad,
    #  fmt=".",
    #  ms=7,
    #  alpha=0.88,
    #  color=colors[1],
    #  elinewidth=1.4,
    #  capsize=1.4,
    #  capthick=1.6,
    #  )
    #  plt.errorbar(
    #  za,
    #  EW_XX_res,
    #  yerr=EW_XX_mad,
    #  fmt=".",
    #  ms=7,
    #  alpha=0.88,
    #  color=colors[2],
    #  elinewidth=1.4,
    #  capsize=1.4,
    #  capthick=1.6,
    #  )
    #  plt.errorbar(
    #  za,
    #  EW_YY_res,
    #  yerr=EW_YY_mad,
    #  fmt=".",
    #  ms=7,
    #  alpha=0.88,
    #  color=colors[3],
    #  elinewidth=1.4,
    #  capsize=1.4,
    #  capthick=1.6,
    #  )

    plt.scatter(
        za,
        NS_XX_res,
        s=16,
        alpha=0.88,
        edgecolor="black",
        linewidth=0.2,
        color=colors[0],
    )
    plt.scatter(
        za,
        NS_YY_res,
        s=16,
        alpha=0.88,
        edgecolor="black",
        linewidth=0.2,
        color=colors[1],
    )
    plt.scatter(
        za,
        EW_XX_res,
        s=16,
        alpha=0.88,
        edgecolor="black",
        linewidth=0.2,
        color=colors[2],
    )
    plt.scatter(
        za,
        EW_YY_res,
        s=16,
        alpha=0.88,
        edgecolor="black",
        linewidth=0.2,
        color=colors[3],
    )

    plt.plot(za, NS_XX_fit, label="NS_XX", linewidth=2, color=colors[0])
    plt.plot(za, NS_YY_fit, label="NS_YY", linewidth=2, color=colors[1])
    plt.plot(za, EW_XX_fit, label="EW_XX", linewidth=2, color=colors[2])
    plt.plot(za, EW_YY_fit, label="EW_YY", linewidth=2, color=colors[3])

    leg = plt.legend(loc="upper right", frameon=True, markerscale=4, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.ylim([-4, 4])
    #  plt.xlabel("Zenith Angle [degrees]")
    plt.ylabel("Residual Power [dB]")
    #  plt.tight_layout()
    #  plt.savefig("ref_rot/rot_resi.pdf")
    return plt


#########################################

if __name__ == "__main__":

    spec, _ = spectral()
    jade, _ = jade()
    #  nside = 512

    ref_model = "ref_rot/ref_dipole_models.npz"

    # Load reference FEE model
    # Rotate the fee models by -pi/2 to move model from spherical (E=0) to Alt/Az (N=0)
    ref_fee_model = np.load(ref_model, allow_pickle=True)

    ref_xx = rotated_fee = rotate_map(
        512, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee_model["XX"]
    )
    ref_yy = rotated_fee = rotate_map(
        512, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee_model["YY"]
    )

    # Plotting stuff
    plt.style.use("seaborn")

    nice_fonts = {
        "font.family": "sans-serif",
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "font.size": 8,
        "legend.fontsize": 5,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }

    plt.rcParams.update(nice_fonts)

    fig = plt.figure(figsize=(3.6, 4.9))

    ax1 = fig.add_subplot(2, 1, 1)
    car_err = slice_residual()
    ax1.set_xticklabels([])
    ax1.set_yticks([-4, -2, 0, 2, 4])
    ax1.set_xlim([-96, 96])
    ax1.set_ylim([-3, 4])
    ax1.text(0.004, 0.93, r"($i$)", horizontalalignment="left", transform=ax1.transAxes)

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.set_yticks([-2, -1, 0, 1, 2])
    ax2.set_xlim([-96, 96])
    ax2.text(0.004, 0.93, r"($ii$)", horizontalalignment="left", transform=ax2.transAxes)
    #  rots = [-30, -10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 30]
    rots = [12, 9, 6, 3, 0]
    #  colors = spec(np.linspace(0.14, 0.63, len(rots)))
    colors = spec([0.14, 0.21, 0.35, 0.63, 0.70])
    legs = [
        Line2D(
            [0], [0], linestyle="None", marker="$NS$", markerfacecolor="k", markersize=9
        ),
        Line2D(
            [0],
            [0],
            linestyle="None",
            marker="$EW$",
            markerfacecolor="k",
            markersize=9.6,
        ),
        Line2D([0], [0], color=colors[0], linewidth=2,),
        Line2D(
            [0],
            [0],
            color=colors[0],
            linewidth=2,
            linestyle="dotted",
            dash_capstyle="round",
        ),
        Line2D([0], [0], color=colors[1], linewidth=2,),
        Line2D(
            [0],
            [0],
            color=colors[1],
            linewidth=2,
            linestyle="dotted",
            dash_capstyle="round",
        ),
        Line2D([0], [0], color=colors[2], linewidth=2,),
        Line2D(
            [0],
            [0],
            color=colors[2],
            linewidth=2,
            linestyle="dotted",
            dash_capstyle="round",
        ),
        Line2D([0], [0], color=colors[3], linewidth=2,),
        Line2D(
            [0],
            [0],
            color=colors[3],
            linewidth=2,
            linestyle="dotted",
            dash_capstyle="round",
        ),
        Line2D([0], [0], color=colors[4], linewidth=2,),
        Line2D(
            [0],
            [0],
            color=colors[4],
            linewidth=2,
            linestyle="dotted",
            dash_capstyle="round",
        ),
    ]
    labels = [
        None,
        None,
        r"$12\degree$",
        r"$12\degree$",
        r"$9\degree$",
        r"$9\degree$",
        r"$6\degree$",
        r"$6\degree$",
        r"$3\degree$",
        r"$3\degree$",
        r"$0\degree$",
        r"$0\degree$",
    ]

    for i, r in enumerate(rots):

        # rotate the ref feko model
        ref_xx_rot = reproject_map(512, r, healpix_array=ref_xx)

        res = ref_xx - ref_xx_rot

        NS, EW = healpix_cardinal_slices(512, res, 85)

        #  plt.plot(NS[1], NS[0], label="NS_XX")
        #  plt.plot(EW[1], EW[0], label="EW_XX")
        plt.plot(
            NS[1][::36],
            poly_fit(NS[1][::36], NS[0][::36], NS[0][::36], 9),
            label=f"NS {r}",
            color=colors[i],
            linewidth=2,
        )

        plt.plot(
            EW[1][::36],
            poly_fit(EW[1][::36], EW[0][::36], EW[0][::36], 9),
            label=f"EW {r}",
            color=colors[i],
            linewidth=2,
            #  dashes=(0.5, 5.),
            linestyle="dotted",
            dash_capstyle="round",
        )

    plt.xlabel("Zenith Angle [degrees]")
    plt.ylabel("Residual Power [dB]")

    ax = fig.add_axes([0.27, 0.314, 0.57, 0.2])
    ax.axis("off")
    leg = ax.legend(legs, labels, mode="expand", ncol=6, frameon=True, loc="upper left")
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.savefig("ref_rot/rot_res_slices_2.pdf")
    plt.close()
