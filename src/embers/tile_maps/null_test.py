"""
Null Tests
----------

A set of tools to perform null tests on reference rf data and reference beam models

"""

from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
from embers.tile_maps.beam_utils import (chisq_fit_gain,
                                         healpix_cardinal_slices, map_slices,
                                         plt_slice, poly_fit, rotate_map)
from matplotlib import pyplot as plt

matplotlib.use("Agg")


def good_ref_maps(nside, map_dir, tile_pair):
    """Collates reference data from 18 good satellites into a good_ref_map

    :param nside: Healpix nside
    :param map_dir: Path to directory with tile_maps_raw, created by :func:`~embers.tile_maps.tile_maps.project_tile_healpix`
    :param tile_pair: List of a pair of mwa tile and reference names. Ex: ["S35XX", "rf0XX"]

    :returns:
        - good_ref_map healpix map, with each pixel of the healpix map containing an array of values from all satellite passes within the pixel

    """

    pointings = ["0", "2", "4", "41"]

    # load data from map .npz file
    f = Path(f"{map_dir}/{tile_pair[0]}_{tile_pair[1]}_sat_maps.npz")
    tile_data = np.load(f, allow_pickle=True)
    tile_data = {key: tile_data[key].item() for key in tile_data}
    ref_map = tile_data["ref_map"]

    # Good sats from which to make plots
    good_sats = [
        25338,
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

    # Empty good ref map
    good_ref_map = [[] for pixel in range(hp.nside2npix(nside))]

    for p in pointings:

        # append to good ref map from all good sat data
        for sat in good_sats:
            for pix in range(hp.nside2npix(nside)):
                good_ref_map[pix].extend(ref_map[p][sat][pix])

    return good_ref_map


def plt_null_test(
    fig=None,
    sub=(None, None, None),
    zen_angle=None,
    del_pow=None,
    del_err=None,
    del_beam=None,
    del_fit=None,
    null_label=None,
    beam_label=None,
    fit_label=None,
    ylabel=True,
    title=None,
):

    """Plot null test between two corresponding slices of reference beam maps

    :param fig: Figure number
    :param sub: Subplot position, tuple of matplotlib indices. Ex: (1, 1, 1)
    :param zen_angle: Array of zenith angles
    :param delta_pow: Null test power b/w corresponding slices of reference beam maps
    :param del_err: Errors on del_pow (MAD)
    :param del_beam: Difference b/w ref beam models, should be zero is both references have the same beam model
    :param del_fit: Polynomial fit to residual power
    :param null_label: Label of null test data
    :param beam_label: Label of beam slice
    :param fit_label: Label of del_fit data
    :param ylabel: If True, plot Y label of plot
    :param title: Plot title

    :returns:
        - ax - :func:`~matplotlib.pyplot.subplot` object

    """

    ax = fig.add_subplot(sub[0], sub[1], sub[2])

    ax.errorbar(
        zen_angle,
        del_pow,
        yerr=del_err,
        fmt=".",
        color="#326765",
        ecolor="#7da87b",
        elinewidth=1.4,
        capsize=1.4,
        capthick=1.6,
        alpha=0.9,
        ms=7,
        label=null_label,
    )

    ax.text(0.02, 0.9, title, horizontalalignment="left", transform=ax.transAxes)

    ax.plot(
        zen_angle, del_beam, color="#c70039", linewidth=1.4, alpha=0.9, label=beam_label
    )

    ax.plot(
        zen_angle, del_fit, color="#ff8264", linewidth=1.4, alpha=0.9, label=fit_label
    )

    ax.set_xlim([-82, 82])
    ax.set_ylim([-10, 10])
    ax.set_xlabel("Zenith Angle [deg]")
    leg = ax.legend(loc="lower center", frameon=True, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    if ylabel is True:
        ax.set_ylabel("Power [dB]")
    else:
        ax.set_yticklabels([])

    return ax


def null_test(nside, za_max, ref_model, map_dir, out_dir):
    """Plot all null tests for reference beam maps

    :param nside: Healpix nside
    :param za_max: Maximum zenith angle
    :param ref_model: Path to feko reference model, saved by :func:`~embers.tile_maps.ref_fee_healpix.ref_healpix_save`
    :param map_dir: Path to directory with tile_maps_raw, created by :func:`~embers.tile_maps.tile_maps.project_tile_healpix`
    :param out_dir: Output directory where null test plots will be saved

    :returns:
        - Null test plot saved to out_dir
        - Reference residuals saved to out_dir

    """

    out_dir.mkdir(parents=True, exist_ok=True)

    tile_pairs = [
        ["S35XX", "rf0XX"],
        ["S35YY", "rf0YY"],
        ["S35XX", "rf1XX"],
        ["S35YY", "rf1YY"],
    ]
    good_rf0XX = rotate_map(
        nside,
        angle=+(1 * np.pi) / 4.0,
        healpix_array=np.asarray(good_ref_maps(nside, map_dir, tile_pairs[0])),
    )
    good_rf0YY = rotate_map(
        nside,
        angle=+(1 * np.pi) / 4.0,
        healpix_array=np.asarray(good_ref_maps(nside, map_dir, tile_pairs[1])),
    )
    good_rf1XX = rotate_map(
        nside,
        angle=+(1 * np.pi) / 4.0,
        healpix_array=np.asarray(good_ref_maps(nside, map_dir, tile_pairs[2])),
    )
    good_rf1YY = rotate_map(
        nside,
        angle=+(1 * np.pi) / 4.0,
        healpix_array=np.asarray(good_ref_maps(nside, map_dir, tile_pairs[3])),
    )

    # NS, EW slices of all four reference tiles
    rf0XX_NS, rf0XX_EW = map_slices(nside, good_rf0XX, za_max)
    rf0YY_NS, rf0YY_EW = map_slices(nside, good_rf0YY, za_max)
    rf1XX_NS, rf1XX_EW = map_slices(nside, good_rf1XX, za_max)
    rf1YY_NS, rf1YY_EW = map_slices(nside, good_rf1YY, za_max)

    # Null test diff in power b/w rf0 & rf1
    ref01_XX_NS = rf0XX_NS[0] - rf1XX_NS[0]
    ref01_XX_EW = rf0XX_EW[0] - rf1XX_EW[0]
    ref01_YY_NS = rf0YY_NS[0] - rf1YY_NS[0]
    ref01_YY_EW = rf0YY_EW[0] - rf1YY_EW[0]

    # Error propogation in null test
    error_ref01_XX_NS = np.sqrt((rf0XX_NS[1]) ** 2 + (rf1XX_NS[1]) ** 2)
    error_ref01_XX_EW = np.sqrt((rf0XX_EW[1]) ** 2 + (rf1XX_EW[1]) ** 2)
    error_ref01_YY_NS = np.sqrt((rf0YY_NS[1]) ** 2 + (rf1YY_NS[1]) ** 2)
    error_ref01_YY_EW = np.sqrt((rf0YY_EW[1]) ** 2 + (rf1YY_EW[1]) ** 2)

    # Load reference FEE model
    ref_fee_model = np.load(ref_model, allow_pickle=True)
    beam_XX = ref_fee_model["XX"]
    beam_YY = ref_fee_model["YY"]

    # Rotate beam models by pi/4 to match rotation of data
    rotated_XX = rotate_map(nside, angle=-(1 * np.pi) / 4.0, healpix_array=beam_XX)
    rotated_YY = rotate_map(nside, angle=-(1 * np.pi) / 4.0, healpix_array=beam_YY)

    # slice the XX rotated map along NS, EW
    XX_NS, XX_EW = healpix_cardinal_slices(nside, rotated_XX, za_max)
    XX_NS_slice, za_NS = XX_NS
    XX_EW_slice, za_EW = XX_EW

    # slice the YY rotated map along NS, EW
    YY_NS, YY_EW = healpix_cardinal_slices(nside, rotated_YY, za_max)
    YY_NS_slice, za_NS = YY_NS
    YY_EW_slice, za_EW = YY_EW

    # Gain offsets for the 8 combinations of data and beam slices
    gain_ref0_XX_NS = chisq_fit_gain(data=rf0XX_NS[0], model=XX_NS_slice)
    gain_ref0_XX_EW = chisq_fit_gain(data=rf0XX_EW[0], model=XX_EW_slice)
    gain_ref1_XX_NS = chisq_fit_gain(data=rf1XX_NS[0], model=XX_NS_slice)
    gain_ref1_XX_EW = chisq_fit_gain(data=rf1XX_EW[0], model=XX_EW_slice)

    gain_ref0_YY_NS = chisq_fit_gain(data=rf0YY_NS[0], model=YY_NS_slice)
    gain_ref0_YY_EW = chisq_fit_gain(data=rf0YY_EW[0], model=YY_EW_slice)
    gain_ref1_YY_NS = chisq_fit_gain(data=rf1YY_NS[0], model=YY_NS_slice)
    gain_ref1_YY_EW = chisq_fit_gain(data=rf1YY_EW[0], model=YY_EW_slice)

    # Scale the data so that it best fits the beam slice
    rf0XX_NS = rf0XX_NS - gain_ref0_XX_NS
    rf0XX_EW = rf0XX_EW - gain_ref0_XX_EW
    rf0YY_NS = rf0YY_NS - gain_ref0_YY_NS
    rf0YY_EW = rf0YY_EW - gain_ref0_YY_EW
    rf1XX_NS = rf1XX_NS - gain_ref1_XX_NS
    rf1XX_EW = rf1XX_EW - gain_ref1_XX_EW
    rf1YY_NS = rf1YY_NS - gain_ref1_YY_NS
    rf1YY_EW = rf1YY_EW - gain_ref1_YY_EW

    # Difference b/w beam model slices.
    # Always 0 because we have only one model for both rf0, rf1
    beam_ref01_XX_NS = XX_NS_slice - XX_NS_slice
    beam_ref01_XX_EW = XX_EW_slice - XX_EW_slice
    beam_ref01_YY_NS = YY_NS_slice - YY_NS_slice
    beam_ref01_YY_EW = YY_EW_slice - YY_EW_slice

    # delta powers
    del_pow_ref0_XX_NS = rf0XX_NS[0] - XX_NS_slice
    del_pow_ref0_XX_EW = rf0XX_EW[0] - XX_EW_slice
    del_pow_ref1_XX_NS = rf1XX_NS[0] - XX_NS_slice
    del_pow_ref1_XX_EW = rf1XX_EW[0] - XX_EW_slice

    del_pow_ref0_YY_NS = rf0YY_NS[0] - YY_NS_slice
    del_pow_ref0_YY_EW = rf0YY_EW[0] - YY_EW_slice
    del_pow_ref1_YY_NS = rf1YY_NS[0] - YY_NS_slice
    del_pow_ref1_YY_EW = rf1YY_EW[0] - YY_EW_slice

    # 3rd order poly fits for residuals
    fit_ref0_XX_NS = poly_fit(za_NS, del_pow_ref0_XX_NS, rf0XX_NS[0], 3)
    fit_ref0_XX_EW = poly_fit(za_EW, del_pow_ref0_XX_EW, rf0XX_EW[0], 3)
    fit_ref1_XX_NS = poly_fit(za_NS, del_pow_ref1_XX_NS, rf1XX_NS[0], 3)
    fit_ref1_XX_EW = poly_fit(za_EW, del_pow_ref1_XX_EW, rf1XX_EW[0], 3)

    fit_ref0_YY_NS = poly_fit(za_NS, del_pow_ref0_YY_NS, rf0YY_NS[0], 3)
    fit_ref0_YY_EW = poly_fit(za_EW, del_pow_ref0_YY_EW, rf0YY_EW[0], 3)
    fit_ref1_YY_NS = poly_fit(za_NS, del_pow_ref1_YY_NS, rf1YY_NS[0], 3)
    fit_ref1_YY_EW = poly_fit(za_EW, del_pow_ref1_YY_EW, rf1YY_EW[0], 3)

    # Difference of power b/w ref0 & ref1 fits
    fit_ref01_XX_NS = fit_ref0_XX_NS - fit_ref1_XX_NS
    fit_ref01_XX_EW = fit_ref0_XX_EW - fit_ref1_XX_EW
    fit_ref01_YY_NS = fit_ref0_YY_NS - fit_ref1_YY_NS
    fit_ref01_YY_EW = fit_ref0_YY_EW - fit_ref1_YY_EW

    ref_res = {
        "za": za_NS,
        "rf0_XX_NS": fit_ref0_XX_NS,
        "rf0_XX_EW": fit_ref0_XX_EW,
        "rf0_YY_NS": fit_ref0_YY_NS,
        "rf0_YY_EW": fit_ref0_YY_EW,
        "rf1_XX_NS": fit_ref1_XX_NS,
        "rf1_XX_EW": fit_ref1_XX_EW,
        "rf1_YY_NS": fit_ref1_YY_NS,
        "rf1_YY_EW": fit_ref1_YY_EW,
    }

    # Save reference residuals to numpy file in out_dir
    # Will be used for errorbars in beam slice plots
    np.save(f"{out_dir}/ref_res", ref_res)

    plt.style.use("seaborn")
    nice_fonts = {
        "font.family": "sans-serif",
        "axes.labelsize": 10,
        "font.size": 10,
        "legend.fontsize": 6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }

    plt.rcParams.update(nice_fonts)

    fig1 = plt.figure(figsize=(12, 7))

    ax1 = plt_slice(
        fig=fig1,
        sub=(3, 4, 1),
        zen_angle=za_NS,
        map_slice=rf0XX_NS[0],
        map_error=rf0XX_NS[1],
        model_slice=XX_NS_slice,
        delta_pow=del_pow_ref0_XX_NS,
        pow_fit=fit_ref0_XX_NS,
        slice_label="ref0XX NS",
        model_label="FEE XX NS",
        title=r"($i$)",
    )

    ax2 = plt_slice(
        fig=fig1,
        sub=(3, 4, 2),
        zen_angle=za_EW,
        map_slice=rf0XX_EW[0],
        map_error=rf0XX_EW[1],
        model_slice=XX_EW_slice,
        delta_pow=del_pow_ref0_XX_EW,
        pow_fit=fit_ref0_XX_EW,
        slice_label="ref0XX EW",
        model_label="FEE XX EW",
        ylabel=False,
        title=r"($ii$)",
    )

    ax3 = plt_slice(
        fig=fig1,
        sub=(3, 4, 5),
        zen_angle=za_NS,
        map_slice=rf1XX_NS[0],
        map_error=rf1XX_NS[1],
        model_slice=XX_NS_slice,
        delta_pow=del_pow_ref1_XX_NS,
        pow_fit=fit_ref1_XX_NS,
        slice_label="ref1XX NS",
        model_label="FEE XX NS",
        title=r"($v$)",
    )

    ax4 = plt_slice(
        fig=fig1,
        sub=(3, 4, 6),
        zen_angle=za_EW,
        map_slice=rf1XX_EW[0],
        map_error=rf1XX_EW[1],
        model_slice=XX_EW_slice,
        delta_pow=del_pow_ref1_XX_EW,
        pow_fit=fit_ref1_XX_EW,
        slice_label="ref1XX EW",
        model_label="FEE XX EW",
        ylabel=False,
        title=r"($vi$)",
    )

    ax5 = plt_null_test(
        fig=fig1,
        sub=(3, 4, 9),
        zen_angle=za_NS,
        del_pow=ref01_XX_NS,
        del_err=error_ref01_XX_NS,
        del_beam=beam_ref01_XX_NS,
        del_fit=fit_ref01_XX_NS,
        null_label="NS rf0-rf1",
        beam_label="FEE Null",
        fit_label="Fit rf0-rf1",
        title=r"($ix$)",
    )

    ax6 = plt_null_test(
        fig=fig1,
        sub=(3, 4, 10),
        zen_angle=za_EW,
        del_pow=ref01_XX_EW,
        del_err=error_ref01_XX_EW,
        del_beam=beam_ref01_XX_EW,
        del_fit=fit_ref01_XX_EW,
        null_label="EW rf0-rf1",
        beam_label="FEE Null",
        fit_label="Fit rf0-rf1",
        ylabel=False,
        title=r"($x$)",
    )

    ax7 = plt_slice(
        fig=fig1,
        sub=(3, 4, 3),
        zen_angle=za_NS,
        map_slice=rf0YY_NS[0],
        map_error=rf0YY_NS[1],
        model_slice=YY_NS_slice,
        delta_pow=del_pow_ref0_YY_NS,
        pow_fit=fit_ref0_YY_NS,
        slice_label="ref0YY NS",
        model_label="FEE YY NS",
        ylabel=False,
        title=r"($iii$)",
    )

    ax8 = plt_slice(
        fig=fig1,
        sub=(3, 4, 4),
        zen_angle=za_EW,
        map_slice=rf0YY_EW[0],
        map_error=rf0YY_EW[1],
        model_slice=YY_EW_slice,
        delta_pow=del_pow_ref0_YY_EW,
        pow_fit=fit_ref0_YY_EW,
        slice_label="ref0YY EW",
        model_label="FEE YY EW",
        ylabel=False,
        title=r"($iv$)",
    )

    ax9 = plt_slice(
        fig=fig1,
        sub=(3, 4, 7),
        zen_angle=za_NS,
        map_slice=rf1YY_NS[0],
        map_error=rf1YY_NS[1],
        model_slice=YY_NS_slice,
        delta_pow=del_pow_ref1_YY_NS,
        pow_fit=fit_ref1_YY_NS,
        slice_label="ref1YY NS",
        model_label="FEE YY NS",
        ylabel=False,
        title=r"($vii$)",
    )

    ax10 = plt_slice(
        fig=fig1,
        sub=(3, 4, 8),
        zen_angle=za_EW,
        map_slice=rf1YY_EW[0],
        map_error=rf1YY_EW[1],
        model_slice=YY_EW_slice,
        delta_pow=del_pow_ref1_YY_EW,
        pow_fit=fit_ref1_YY_EW,
        slice_label="ref1YY EW",
        model_label="FEE YY EW",
        ylabel=False,
        title=r"($viii$)",
    )

    ax11 = plt_null_test(
        fig=fig1,
        sub=(3, 4, 11),
        zen_angle=za_NS,
        del_pow=ref01_YY_NS,
        del_err=error_ref01_YY_NS,
        del_beam=beam_ref01_YY_NS,
        del_fit=fit_ref01_YY_NS,
        null_label="NS rf0-rf1",
        beam_label="FEE Null",
        fit_label="Fit rf0-rf1",
        ylabel=False,
        title=r"($xi$)",
    )

    ax12 = plt_null_test(
        fig=fig1,
        sub=(3, 4, 12),
        zen_angle=za_EW,
        del_pow=ref01_YY_EW,
        del_err=error_ref01_YY_EW,
        del_beam=beam_ref01_YY_EW,
        del_fit=fit_ref01_YY_EW,
        null_label="EW rf0-rf1",
        beam_label="FEE Null",
        fit_label="Fit rf0-rf1",
        ylabel=False,
        title=r"($xii$)",
    )

    plt.tight_layout()
    fig1.savefig(f"{out_dir}/null_test.pdf", bbox_inches="tight")
