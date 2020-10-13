"""
NULL SATS.
----------
"""

import argparse
from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
from embers.rf_tools.colormaps import spectral
from matplotlib import pyplot as plt
from numpy.polynomial import polynomial as poly
from scipy.stats import median_absolute_deviation as mad

matplotlib.use("Agg")
_spec, _ = spectral()
colors = _spec([0.72, 0.30, 0.18])

parser = argparse.ArgumentParser(
    description="""
    Null Sats paper plot
    """
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots",
    help="Output directory. Default=../embers_out/paper_plots",
)
parser.add_argument(
    "--map_dir",
    metavar="\b",
    default="../embers_out/tile_maps/tile_maps/tile_maps_raw/",
    help="Output directory. Default=../embers_out/tile_maps/tile_maps/tile_maps_raw/",
)
parser.add_argument(
    "--ref_model",
    metavar="\b",
    default="../embers_out/tile_maps/ref_models/ref_dipole_models.npz",
    help="Healpix reference FEE model file. default=../embers_out/tile_maps/ref_models/ref_dipole_models.npz",
)
parser.add_argument(
    "--nside", metavar="\b", type=int, default=32, help="Healpix Nside. Default = 32",
)

args = parser.parse_args()

out_dir = Path(args.out_dir)
map_dir = Path(args.map_dir)
ref_model = args.ref_model
nside = args.nside

out_dir.mkdir(parents=True, exist_ok=True)


def ring_indices(nside=None):
    """Healpix pix indices of NS, EW slices and above horizon"""

    # theta phi values of each pixel
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    za = np.arange(0, 87, 2)

    ring_indices = []
    for i in range(len(za)):
        if i < 43:
            r_indices = np.where(
                np.logical_and(θ >= np.radians(za[i]), θ < np.radians(za[i + 1]))
            )[0]
            ring_indices.append(r_indices)
    za_cen = za[:-1] + 1
    return [za_cen, ring_indices]


def good_maps(ref_map):
    """Creates a ref map with only good satellites"""

    pointings = ["0", "2", "4"]

    # load data from map .npz file
    f = Path(f"{map_dir}/{ref_map}")
    tile_data = np.load(f, allow_pickle=True)
    tile_data = {key: tile_data[key].item() for key in tile_data}
    ref_map = tile_data["ref_map"]

    # Good sats from which to make plots
    good_sats = [
        25338,
        25982,
        25984,
        25985,
        28654,
        40086,
        40087,
        40091,
        41179,
        41180,
        41182,
        41183,
        41184,
        41185,
        41187,
        41188,
        41189,
        44387,
    ]

    orbcomm = [
        25982,
        25984,
        25985,
        40086,
        40087,
        40091,
        41179,
        41180,
        41182,
        41183,
        41184,
        41185,
        41187,
        41188,
        41189,
    ]
    noaa = [25338, 28654]

    meteor = [44387]

    sat_types = [orbcomm, noaa, meteor]
    sat_types_names = ["orbcomm", "noaa", "meteor"]

    map_dict = {}

    for index, s_type in enumerate(sat_types):

        # Empty good map
        good_map = [[] for pixel in range(hp.nside2npix(nside))]

        for p in pointings:

            # append to good map from all good sat data
            for sat in s_type:
                for pix in range(hp.nside2npix(nside)):
                    good_map[pix].extend(ref_map[p][sat][pix])

        mad_map = []
        for j in good_map:
            if j != []:
                j = np.asarray(j)
                j = j[~np.isnan(j)]
                mad_map.append(mad(j))
            else:
                mad_map.append(np.nan)

        good_map = [np.nanmedian(pixel) for pixel in good_map]

        map_type = [good_map, mad_map]

        map_dict[sat_types_names[index]] = map_type

    return map_dict


def poly_fit(x, y, order):
    """Fit polynominal of order to data"""

    x = np.asarray(x)
    y = np.asarray(y)

    coefs = poly.polyfit(x, y, order)
    fit = poly.polyval(x, coefs)
    return fit


try:

    ref_tiles = [
        "S35XX_rf0XX_sat_maps.npz",
        "S35YY_rf0YY_sat_maps.npz",
        "S35XX_rf1XX_sat_maps.npz",
        "S35YY_rf1YY_sat_maps.npz",
    ]

    za, indices = ring_indices(nside=nside)

    ref_prof = []

    sat_types_names = ["orbcomm", "noaa", "meteor"]
    for s_type in sat_types_names:

        # median and mad values of induvidual ref maps
        good_rf0XX, mad_rf0XX = good_maps(ref_tiles[0])[s_type]
        good_rf0YY, mad_rf0YY = good_maps(ref_tiles[1])[s_type]
        good_rf1XX, mad_rf1XX = good_maps(ref_tiles[2])[s_type]
        good_rf1YY, mad_rf1YY = good_maps(ref_tiles[3])[s_type]

        # Load reference FEE model
        ref_fee_model = np.load(ref_model, allow_pickle=True)
        beam_XX = ref_fee_model["XX"]
        beam_YY = ref_fee_model["YY"]

        # residuals of individual maps
        res_rf0XX = good_rf0XX - beam_XX
        res_rf0YY = good_rf0YY - beam_YY
        res_rf1XX = good_rf1XX - beam_XX
        res_rf1YY = good_rf1YY - beam_YY

        sum_array = np.array([res_rf0XX, res_rf0YY, res_rf1XX, res_rf1YY])
        mad_array = np.array([mad_rf0XX, mad_rf0YY, mad_rf1XX, mad_rf1YY])

        # sum along corresponding pixels
        res_median = np.nanmedian(sum_array, axis=0)
        mad_average = np.nanmean(mad_array, axis=0)

        # median of values in zenith angle rings
        res_rings = []
        mad_rings = []
        for i in indices:
            res_rings.append(np.nanmedian(res_median[i]))
            mad_rings.append(np.nanmedian(mad_average[i]))

        # Scale the residuals to 0 median
        res_med = np.median(res_rings)
        res_rings = np.array(res_rings) - res_med
        fit_res = poly_fit(za, res_rings, 8)

        prof_data = [res_rings, mad_rings, fit_res]

        ref_prof.append(prof_data)

    plt.style.use("seaborn")
    nice_fonts = {
        # Use LaTeX to write all text
        # "text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 8,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }

    plt.rcParams.update(nice_fonts)

    fig = plt.figure(figsize=(3.6, 2.4))

    #  plt.style.use("seaborn")
    # plt.errorbar(
    #       za, ref_prof[0][0], yerr=ref_prof[0][1],
    #       fmt='.', color='#326765', ecolor='#7da87b',
    #       elinewidth=1.4, capsize=1.4, capthick=1.6,
    #       alpha=0.9, ms=9, label=r'$\Delta$ref [$\theta$]')
    #
    # plt.plot(za, ref_prof[0][2], linewidth=1.4, alpha=1, color='#fa4659', label=r'8$^{th}$ order fit')
    #
    # plt.errorbar(
    #       za, ref_prof[1][0], yerr=ref_prof[1][1],
    #       fmt='.', color='#326765', ecolor='#7da87b',
    #       elinewidth=1.4, capsize=1.4, capthick=1.6,
    #       alpha=0.9, ms=9, label=r'$\Delta$ref [$\theta$]')
    #
    # plt.plot(za, ref_prof[1][2], linewidth=1.4, alpha=1, color='#fa4659', label=r'8$^{th}$ order fit')
    #
    # plt.errorbar(
    #       za, ref_prof[2][0], yerr=ref_prof[2][1],
    #       fmt='.', color='#326765', ecolor='#7da87b',
    #       elinewidth=1.4, capsize=1.4, capthick=1.6,
    #       alpha=0.9, ms=9, label=r'$\Delta$ref [$\theta$]')
    #
    # plt.plot(za, ref_prof[2][2], linewidth=1.4, alpha=1, color='#fa4659', label=r'8$^{th}$ order fit')

    #  plt.plot(
    #  za, ref_prof[0][2], linewidth=1.4, alpha=1, color="#070d59", label="ORBCOMM"
    #  )
    #  plt.scatter(za, ref_prof[0][0], color="#070d59", alpha=0.6, marker="o", s=14)

    #  plt.plot(
    #  za, ref_prof[2][2], linewidth=1.4, alpha=1, color="#5893d4", label="METEOR"
    #  )
    #  plt.scatter(za, ref_prof[2][0], color="#5893d4", alpha=0.6, marker="o", s=14)

    #  plt.plot(za, ref_prof[1][2], linewidth=1.4, alpha=1, color="#729d39", label="NOAA")
    #  plt.scatter(za, ref_prof[1][0], color="#729d39", alpha=0.6, marker="o", s=14)

    plt.plot(
        za, ref_prof[0][2], linewidth=1.4, alpha=1, color=colors[0], label="ORBCOMM"
    )
    plt.scatter(za, ref_prof[0][0], color=colors[0], alpha=0.6, marker="o", s=14)

    plt.plot(
        za, ref_prof[2][2], linewidth=1.4, alpha=1, color=colors[1], label="METEOR"
    )
    plt.scatter(za, ref_prof[2][0], color=colors[1], alpha=0.6, marker="o", s=14)

    plt.plot(za, ref_prof[1][2], linewidth=1.4, alpha=1, color=colors[2], label="NOAA")
    plt.scatter(za, ref_prof[1][0], color=colors[2], alpha=0.6, marker="o", s=14)

    leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.xlabel(r"Zenith Angle [degrees]")
    plt.ylabel("Residual Power [dB]")
    plt.tight_layout()

    plt.savefig(f"{out_dir}/ref_residuals.pdf", bbox_inches="tight")

    print(f"NULL SATS saved to {out_dir}")


except Exception as e:
    print(e)
