from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
from embers.rf_tools.colormaps import spectral
from embers.tile_maps.beam_utils import (chisq_fit_gain,
                                         healpix_cardinal_slices, map_slices,
                                         poly_fit, rotate_map)
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.use("Agg")
cmap, _ = spectral()


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


def plt_slice(
    fig=None,
    sub=None,
    zen_angle=None,
    map_slice=None,
    map_error=None,
    model_slice=None,
    delta_pow=None,
    pow_fit=None,
    slice_label=None,
    model_label=None,
    ylabel=True,
):

    """Plot a slice of the beam, with measured
    data, errorbars, and fit the simulated beam
    to the data. Also plot the diff b/w data and
    the model"""

    ax = fig.add_subplot(sub)

    ax.errorbar(
        zen_angle,
        map_slice,
        yerr=map_error,
        fmt=".",
        color="#326765",
        ecolor="#7da87b",
        elinewidth=1.2,
        capsize=1.2,
        capthick=1.4,
        alpha=0.9,
        label=slice_label,
    )

    ax.plot(
        zen_angle,
        model_slice,
        color="#c70039",
        linewidth=1.2,
        alpha=0.9,
        label=model_label,
    )

    # ax.set_ylim(bottom=-30)
    # ax.set_xlabel('Zenith Angle (degrees)')
    # ax.legend(loc='lower center')
    leg = ax.legend(loc="lower center", frameon=True, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)
    ax.set_xlim([-82, 82])
    ax.set_ylim([-26, 12])

    divider = make_axes_locatable(ax)
    dax = divider.append_axes("bottom", size="30%", pad=0.1)

    # dax = fig.add_subplot(2,1,2)
    dax.scatter(zen_angle, delta_pow, marker=".", color="#27296d")
    dax.plot(zen_angle, pow_fit, linewidth=1.2, alpha=0.9, color="#ff8264")
    dax.set_xlim([-82, 82])
    dax.set_xticklabels([])
    dax.set_ylim([-5, 5])
    if ylabel is True:
        ax.set_ylabel("Power [dB]")
        dax.set_ylabel(r"$\Delta ref$ [dB]")
        # dax.set_ylabel('Data - Model (dB)')
    else:
        dax.set_yticklabels([])
        ax.set_yticklabels([])

    return ax


def plt_null_test(
    fig=None,
    sub=None,
    zen_angle=None,
    del_pow=None,
    del_err=None,
    del_beam=None,
    del_fit=None,
    null_label=None,
    beam_label=None,
    fit_label=None,
    ylabel=True,
):

    """Plot graphs for a null test"""

    ax = fig.add_subplot(sub)

    ax.errorbar(
        zen_angle,
        del_pow,
        yerr=del_err,
        fmt=".",
        color="#326765",
        ecolor="#7da87b",
        elinewidth=1.2,
        capsize=1.2,
        capthick=1.4,
        alpha=0.9,
        label=null_label,
    )

    ax.plot(
        zen_angle, del_beam, color="#c70039", linewidth=1.2, alpha=0.9, label=beam_label
    )

    ax.plot(
        zen_angle, del_fit, color="#ff8264", linewidth=1.2, alpha=0.9, label=fit_label
    )

    ax.set_xlim([-82, 82])
    ax.set_ylim([-10, 10])
    ax.set_xlabel("Zenith Angle [degrees]")
    # ax.legend(loc='upper left')
    leg = ax.legend(loc="upper left", frameon=True, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    if ylabel is True:
        ax.set_ylabel("Power [dB]")
    else:
        ax.set_yticklabels([])

    return ax


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Plot healpix map of reference data
        """
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="../../outputs/tile_maps/null_test/",
        help="Output directory. Default=../../outputs/tile_maps/null_test/",
    )

    parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="../../outputs/tile_maps/tile_maps_raw/",
        help="Output directory. Default=../../outputs/tile_maps/tile_maps_raw/",
    )

    parser.add_argument(
        "--ref_model",
        metavar="\b",
        default="../../outputs/reproject_ref/ref_dipole_models.npz",
        help="Healpix reference FEE model file. default=../../outputs/reproject_ref/ref_dipole_models.npz",
    )

    parser.add_argument(
        "--nside",
        metavar="\b",
        type=int,
        default=32,
        help="Healpix Nside. Default = 32",
    )

    parser.add_argument(
        "--za_max",
        metavar="\b",
        type=int,
        default=80,
        help="Maximum zenith angle. Default = 80 deg",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    map_dir = Path(args.map_dir)
    ref_model = args.ref_model
    nside = args.nside
    za_max = args.za_max

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

    plt.style.use("seaborn")
    fig1 = plt.figure(figsize=(8, 10))

    ax1 = plt_slice(
        fig=fig1,
        sub=321,
        zen_angle=za_NS,
        map_slice=rf0XX_NS[0],
        map_error=rf0XX_NS[1],
        model_slice=XX_NS_slice,
        delta_pow=del_pow_ref0_XX_NS,
        pow_fit=fit_ref0_XX_NS,
        slice_label="ref0XX NS",
        model_label="FEE XX NS",
    )

    ax2 = plt_slice(
        fig=fig1,
        sub=322,
        zen_angle=za_EW,
        map_slice=rf0XX_EW[0],
        map_error=rf0XX_EW[1],
        model_slice=XX_EW_slice,
        delta_pow=del_pow_ref0_XX_EW,
        pow_fit=fit_ref0_XX_EW,
        slice_label="ref0XX EW",
        model_label="FEE XX EW",
        ylabel=False,
    )

    ax3 = plt_slice(
        fig=fig1,
        sub=323,
        zen_angle=za_NS,
        map_slice=rf1XX_NS[0],
        map_error=rf1XX_NS[1],
        model_slice=XX_NS_slice,
        delta_pow=del_pow_ref1_XX_NS,
        pow_fit=fit_ref1_XX_NS,
        slice_label="ref1XX NS",
        model_label="FEE XX NS",
    )

    ax4 = plt_slice(
        fig=fig1,
        sub=324,
        zen_angle=za_EW,
        map_slice=rf1XX_EW[0],
        map_error=rf1XX_EW[1],
        model_slice=XX_EW_slice,
        delta_pow=del_pow_ref1_XX_EW,
        pow_fit=fit_ref1_XX_EW,
        slice_label="ref1XX EW",
        model_label="FEE XX EW",
        ylabel=False,
    )

    ax5 = plt_null_test(
        fig=fig1,
        sub=325,
        zen_angle=za_NS,
        del_pow=ref01_XX_NS,
        del_err=error_ref01_XX_NS,
        del_beam=beam_ref01_XX_NS,
        del_fit=fit_ref01_XX_NS,
        null_label="NS rf0-rf1",
        beam_label="FEE Null",
        fit_label="Fit rf0-rf1",
    )

    ax6 = plt_null_test(
        fig=fig1,
        sub=326,
        zen_angle=za_EW,
        del_pow=ref01_XX_EW,
        del_err=error_ref01_XX_EW,
        del_beam=beam_ref01_XX_EW,
        del_fit=fit_ref01_XX_EW,
        null_label="EW rf0-rf1",
        beam_label="FEE Null",
        fit_label="Fit rf0-rf1",
        ylabel=False,
    )

    plt.tight_layout()
    fig1.savefig(f"{out_dir}/null_test_XX_slices.png")

    fig2 = plt.figure(figsize=(8, 10))

    ax7 = plt_slice(
        fig=fig2,
        sub=321,
        zen_angle=za_NS,
        map_slice=rf0YY_NS[0],
        map_error=rf0YY_NS[1],
        model_slice=YY_NS_slice,
        delta_pow=del_pow_ref0_YY_NS,
        pow_fit=fit_ref0_YY_NS,
        slice_label="ref0YY NS",
        model_label="FEE YY NS",
    )

    ax8 = plt_slice(
        fig=fig2,
        sub=322,
        zen_angle=za_EW,
        map_slice=rf0YY_EW[0],
        map_error=rf0YY_EW[1],
        model_slice=YY_EW_slice,
        delta_pow=del_pow_ref0_YY_EW,
        pow_fit=fit_ref0_YY_EW,
        slice_label="ref0YY EW",
        model_label="FEE YY EW",
        ylabel=False,
    )

    ax9 = plt_slice(
        fig=fig2,
        sub=323,
        zen_angle=za_NS,
        map_slice=rf1YY_NS[0],
        map_error=rf1YY_NS[1],
        model_slice=YY_NS_slice,
        delta_pow=del_pow_ref1_YY_NS,
        pow_fit=fit_ref1_YY_NS,
        slice_label="ref1YY NS",
        model_label="FEE YY NS",
    )

    ax10 = plt_slice(
        fig=fig2,
        sub=324,
        zen_angle=za_EW,
        map_slice=rf1YY_EW[0],
        map_error=rf1YY_EW[1],
        model_slice=YY_EW_slice,
        delta_pow=del_pow_ref1_YY_EW,
        pow_fit=fit_ref1_YY_EW,
        slice_label="ref1YY EW",
        model_label="FEE YY EW",
        ylabel=False,
    )

    ax11 = plt_null_test(
        fig=fig2,
        sub=325,
        zen_angle=za_NS,
        del_pow=ref01_YY_NS,
        del_err=error_ref01_YY_NS,
        del_beam=beam_ref01_YY_NS,
        del_fit=fit_ref01_YY_NS,
        null_label="NS rf0-rf1",
        beam_label="FEE Null",
        fit_label="Fit rf0-rf1",
    )

    ax12 = plt_null_test(
        fig=fig2,
        sub=326,
        zen_angle=za_EW,
        del_pow=ref01_YY_EW,
        del_err=error_ref01_YY_EW,
        del_beam=beam_ref01_YY_EW,
        del_fit=fit_ref01_YY_EW,
        null_label="EW rf0-rf1",
        beam_label="FEE Null",
        fit_label="Fit rf0-rf1",
        ylabel=False,
    )

    plt.tight_layout()
    fig2.savefig(f"{out_dir}/null_test_YY_slices.png")
