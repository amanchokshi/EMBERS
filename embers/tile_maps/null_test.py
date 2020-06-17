import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly


def hp_slices_horizon(nside=None):
    """Healpix pix indices of NS, EW slices and above horizon"""
    # theta phi values of each pixel
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    # healpix indices above the horizon
    above_horizon_indices = np.where(θ <= np.radians(80))[0]

    # pixel coords above the horizon
    ɸ_above_horizon = ɸ[above_horizon_indices]

    NS_indices = []
    EW_indices = []

    # pixel indices along N, E, S, W slices
    # order the indices such that they proceed from N -> S or E -> W
    n_slice = sorted(
        np.where((np.round(np.degrees(ɸ_above_horizon))) == 45)[0], reverse=True
    )
    e_slice = sorted(
        np.where((np.round(np.degrees(ɸ_above_horizon))) == 135)[0], reverse=True
    )
    s_slice = sorted(np.where((np.round(np.degrees(ɸ_above_horizon))) == 225)[0])
    w_slice = sorted(np.where((np.round(np.degrees(ɸ_above_horizon))) == 315)[0])

    NS_indices.extend(n_slice)
    NS_indices.extend(s_slice)
    EW_indices.extend(e_slice)
    EW_indices.extend(w_slice)

    return [NS_indices, EW_indices, above_horizon_indices]


def slice_map(hp_map):
    """slices healpix map along NS, EW"""

    NS_indices, EW_indices, _ = hp_slices_horizon(nside)

    θ_NS, ɸ_NS = np.degrees(hp.pix2ang(nside, NS_indices))
    θ_EW, ɸ_EW = np.degrees(hp.pix2ang(nside, EW_indices))

    zenith_angle_NS = []
    for i, j in zip(θ_NS, ɸ_NS):
        if j <= 180:
            zenith_angle_NS.append(-1 * i)
        else:
            zenith_angle_NS.append(i)

    zenith_angle_EW = []
    for i, j in zip(θ_EW, ɸ_EW):
        if j <= 180:
            zenith_angle_EW.append(-1 * i)
        else:
            zenith_angle_EW.append(i)

    NS_data = [hp_map[NS_indices], zenith_angle_NS]
    EW_data = [hp_map[EW_indices], zenith_angle_EW]

    return [NS_data, EW_data]


def nan_mad(ref_map):
    """Compute mad while ignoring nans"""
    ref_map_mad = []
    for j in ref_map:
        if j != []:
            j = np.asarray(j)
            j = j[~np.isnan(j)]
            ref_map_mad.append(mad(j))
        else:
            ref_map_mad.append(np.nan)

    ref_map_mad = np.asarray(ref_map_mad)
    ref_map_mad[np.where(ref_map_mad == np.nan)] = np.nanmean(ref_map_mad)

    return ref_map_mad


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

    # Empty good map
    good_map = [[] for pixel in range(hp.nside2npix(nside))]

    for p in pointings:

        # append to good map from all good sat data
        for sat in good_sats:
            for pix in range(hp.nside2npix(nside)):
                good_map[pix].extend(ref_map[p][sat][pix])

    return good_map


def ref_map_slice(good_map):
    """slices ref healpix map along NS & EW"""

    ref_map_NS, ref_map_EW = slice_map(np.asarray(good_map))

    ref_med_map_NS = np.asarray(
        [(np.nanmedian(i) if i != [] else np.nan) for i in ref_map_NS[0]]
    )
    # Scale mean map such that the max value is 0
    ref_med_map_scaled_NS = np.asarray(
        [i - np.nanmax(ref_med_map_NS) for i in ref_med_map_NS]
    )
    # ref_mad_map_NS = np.asarray([mad(i) for i in ref_map_NS[0]])
    ref_mad_map_NS = np.asarray(nan_mad(ref_map_NS[0]))
    za_NS = ref_map_NS[1]

    ref_med_map_EW = np.asarray(
        [(np.nanmedian(i) if i != [] else np.nan) for i in ref_map_EW[0]]
    )
    # Scale mean map such that the max value is 0
    ref_med_map_scaled_EW = np.asarray(
        [i - np.nanmax(ref_med_map_EW) for i in ref_med_map_EW]
    )
    # ref_mad_map_EW = np.asarray([mad(i) for i in ref_map_EW[0]])
    ref_mad_map_EW = np.asarray(nan_mad(ref_map_EW[0]))
    za_EW = ref_map_EW[1]

    NS_data = [ref_med_map_scaled_NS, ref_mad_map_NS, za_NS]
    EW_data = [ref_med_map_scaled_EW, ref_mad_map_EW, za_EW]

    return [NS_data, EW_data]


# rotate func written by Jack Line
def rotate(nside, angle=None, healpix_array=None, savetag=None, flip=False):
    """Takes in a healpix array, rotates it by the desired angle, and saves it.
    Optionally flip the data, changes east-west into west-east because
    astronomy"""

    # theta phi values of each pixel
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    new_hp_inds = hp.ang2pix(nside, θ, ɸ + angle)

    ##Flip the data to match astro conventions
    if flip == True:
        new_angles = []
        for phi in ɸ:
            if phi <= np.pi:
                new_angles.append(np.pi - phi)
            else:
                new_angles.append(3 * np.pi - phi)
        new_hp_inds = hp.ang2pix(nside, ɸ, np.asarray(new_angles))

    ##Save the array in the new order
    if savetag:
        np.savez_compressed(savetag, beammap=healpix_array[new_hp_inds])

    return healpix_array[new_hp_inds]


# chisquared minimization to best fit map to data
def fit_gain(map_data=None, map_error=None, beam=None):
    """Fit the beam model to the measured data using
    chisquared minimization"""

    bad_values = np.isnan(map_data)
    map_data = map_data[~bad_values]
    map_error = map_error[~bad_values]

    map_error[np.where(map_error == 0)] = np.mean(map_error)

    def chisqfunc(gain):
        model = beam[~bad_values] + gain
        chisq = sum((map_data - model) ** 2)
        # chisq = sum(((map_data - model)/map_error)**2)
        return chisq

    x0 = np.array([0])

    result = opt.minimize(chisqfunc, x0)

    return result.x


def poly_fit(x, y, map_data, order):
    """Fit polynominal of order to data"""

    x = np.asarray(x)
    y = np.asarray(y)

    bad_values = np.isnan(map_data)
    x_good = x[~bad_values]
    y_good = y[~bad_values]
    coefs = poly.polyfit(x_good, y_good, order)
    fit = poly.polyval(x, coefs)
    return fit


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
    for l in leg.legendHandles:
        l.set_alpha(1)
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
    for l in leg.legendHandles:
        l.set_alpha(1)

    if ylabel is True:
        ax.set_ylabel("Power [dB]")
    else:
        ax.set_yticklabels([])

    return ax


if __name__ == "__main__":

    import argparse
    import numpy as np
    from pathlib import Path
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gs
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.stats import median_absolute_deviation as mad

    import sys

    sys.path.append("../decode_rf_data")
    from plot_tile_maps import plot_healpix
    from colormap import spectral

    # Custom spectral colormap
    cmap = spectral()

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

    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    map_dir = Path(args.map_dir)
    ref_model = args.ref_model
    nside = args.nside

    out_dir.mkdir(parents=True, exist_ok=True)

    ref_tiles = [
        "S35XX_rf0XX_sat_maps.npz",
        "S35YY_rf0YY_sat_maps.npz",
        "S35XX_rf1XX_sat_maps.npz",
        "S35YY_rf1YY_sat_maps.npz",
    ]
    good_rf0XX = rotate(
        nside,
        angle=+(1 * np.pi) / 4.0,
        healpix_array=np.asarray(good_maps(ref_tiles[0])),
    )
    good_rf0YY = rotate(
        nside,
        angle=+(1 * np.pi) / 4.0,
        healpix_array=np.asarray(good_maps(ref_tiles[1])),
    )
    good_rf1XX = rotate(
        nside,
        angle=+(1 * np.pi) / 4.0,
        healpix_array=np.asarray(good_maps(ref_tiles[2])),
    )
    good_rf1YY = rotate(
        nside,
        angle=+(1 * np.pi) / 4.0,
        healpix_array=np.asarray(good_maps(ref_tiles[3])),
    )

    # NS, EW slices of all four reference tiles
    rf0XX_NS, rf0XX_EW = ref_map_slice(good_rf0XX)
    rf0YY_NS, rf0YY_EW = ref_map_slice(good_rf0YY)
    rf1XX_NS, rf1XX_EW = ref_map_slice(good_rf1XX)
    rf1YY_NS, rf1YY_EW = ref_map_slice(good_rf1YY)

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
    rotated_XX = rotate(nside, angle=-(1 * np.pi) / 4.0, healpix_array=beam_XX)
    rotated_YY = rotate(nside, angle=-(1 * np.pi) / 4.0, healpix_array=beam_YY)

    # These plots show that the pi/4 rotation was correct
    # plot_healpix(data_map=rotated_XX,sub=(1,1,1), cmap=cmap, vmin=-40, vmax=-20)
    # plot_healpix(data_map=rotated_YY,sub=(1,1,1), cmap=cmap, vmin=-40, vmax=-20)
    # plt.show()

    # slice the XX rotated map along NS, EW
    XX_NS, XX_EW = slice_map(rotated_XX)
    XX_NS_slice, za_NS = XX_NS
    XX_EW_slice, za_EW = XX_EW

    # slice the YY rotated map along NS, EW
    YY_NS, YY_EW = slice_map(rotated_YY)
    YY_NS_slice, za_NS = YY_NS
    YY_EW_slice, za_EW = YY_EW

    # Gain offsets for the 8 combinations of data and beam slices
    gain_ref0_XX_NS = fit_gain(
        map_data=rf0XX_NS[0], map_error=rf0XX_NS[1], beam=XX_NS_slice
    )
    gain_ref0_XX_EW = fit_gain(
        map_data=rf0XX_EW[0], map_error=rf0XX_EW[1], beam=XX_EW_slice
    )
    gain_ref1_XX_NS = fit_gain(
        map_data=rf1XX_NS[0], map_error=rf1XX_NS[1], beam=XX_NS_slice
    )
    gain_ref1_XX_EW = fit_gain(
        map_data=rf1XX_EW[0], map_error=rf1XX_EW[1], beam=XX_EW_slice
    )

    gain_ref0_YY_NS = fit_gain(
        map_data=rf0YY_NS[0], map_error=rf0YY_NS[1], beam=YY_NS_slice
    )
    gain_ref0_YY_EW = fit_gain(
        map_data=rf0YY_EW[0], map_error=rf0YY_EW[1], beam=YY_EW_slice
    )
    gain_ref1_YY_NS = fit_gain(
        map_data=rf1YY_NS[0], map_error=rf1YY_NS[1], beam=YY_NS_slice
    )
    gain_ref1_YY_EW = fit_gain(
        map_data=rf1YY_EW[0], map_error=rf1YY_EW[1], beam=YY_EW_slice
    )

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
