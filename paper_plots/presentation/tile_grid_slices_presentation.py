import sys
import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly
from mwa_pb.mwa_sweet_spots import all_grid_points

sys.path.append("../../decode_rf_data")
from colormap import spectral, jade, kelp

# Custom spectral colormap
jade, _ = jade()


def hp_slices_horizon(nside=None, zenith_angle=90):
    """Healpix pix indices of NS, EW slices and above horizon"""
    # theta phi values of each pixel
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    # healpix indices above the horizon
    above_horizon_indices = np.where(θ <= np.radians(zenith_angle))[0]

    # pixel coords above the horizon
    ɸ_above_horizon = ɸ[above_horizon_indices]

    NS_indices = []
    EW_indices = []

    # pixel indices along N, E, S, W slices
    # order the indices such that they proceed from N -> S or E -> W
    n_slice = sorted(np.where((np.round(np.degrees(ɸ_above_horizon))) == 45)[0])
    e_slice = sorted(np.where((np.round(np.degrees(ɸ_above_horizon))) == 135)[0])
    s_slice = sorted(
        np.where((np.round(np.degrees(ɸ_above_horizon))) == 225)[0], reverse=True
    )
    w_slice = sorted(
        np.where((np.round(np.degrees(ɸ_above_horizon))) == 315)[0], reverse=True
    )

    NS_indices.extend(s_slice)
    NS_indices.extend(n_slice)
    EW_indices.extend(w_slice)
    EW_indices.extend(e_slice)

    return [NS_indices, EW_indices, above_horizon_indices]


def slice_map(hp_map, za):
    """slices healpix map along NS, EW"""

    NS_indices, EW_indices, _ = hp_slices_horizon(nside, zenith_angle=za)

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


def tile_map_slice(good_map, za):
    """slices ref healpix map along NS & EW"""

    ref_map_NS, ref_map_EW = slice_map(np.asarray(good_map), za)

    ref_med_map_NS = np.asarray(
        [(np.nanmedian(i) if i != [] else np.nan) for i in ref_map_NS[0]]
    )
    ref_mad_map_NS = np.asarray(nan_mad(ref_map_NS[0]))
    za_NS = ref_map_NS[1]

    ref_med_map_EW = np.asarray(
        [(np.nanmedian(i) if i != [] else np.nan) for i in ref_map_EW[0]]
    )
    ref_mad_map_EW = np.asarray(nan_mad(ref_map_EW[0]))
    za_EW = ref_map_EW[1]

    NS_data = [ref_med_map_NS, ref_mad_map_NS, za_NS]
    EW_data = [ref_med_map_EW, ref_mad_map_EW, za_EW]

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
    sub=[1, 1, 1],
    zen_angle=None,
    map_slice=None,
    map_error=None,
    model_slice=None,
    delta_pow=None,
    pow_fit=None,
    slice_label=None,
    model_label=None,
    xlabel=False,
    ylabel=False,
    title=None,
):

    """Plot a slice of the beam, with measured
    data, errorbars, and fit the simulated beam
    to the data. Also plot the diff b/w data and
    the model"""

    ax = fig.add_subplot(sub[0], sub[1], sub[2])

    ax.errorbar(
        zen_angle,
        map_slice,
        yerr=map_error,
        fmt=".",
        color="#a8df65",
        ecolor="#7da87b",
        elinewidth=2.1,
        capsize=2.1,
        capthick=2.6,
        alpha=0.9,
        label=slice_label,
    )

    ax.plot(
        zen_angle,
        model_slice,
        color="#c70039",
        linewidth=2.1,
        alpha=0.9,
        label=model_label,
    )
    ax.set_xticks([-75, -50, -25, 0, 25, 50, 75])

    leg = ax.legend(loc="lower center", frameon=True, framealpha=0.3, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for l in leg.legendHandles:
        l.set_alpha(1)
    for text in leg.get_texts():
        plt.setp(text, color="w")

    # ax.set_ylim(bottom=-30)
    ax.set_xlim([-92, 92])
    ax.set_ylim([-50, 5])
    ax.set_title(title, loc="left", color="w")

    divider = make_axes_locatable(ax)
    dax = divider.append_axes("bottom", size="30%", pad=0.05)
    dax.set_xticks([-75, -50, -25, 0, 25, 50, 75])

    if xlabel is True:
        dax.set_xlabel("Zenith Angle [deg]")
    else:
        dax.set_xticklabels([])
        ax.set_xticklabels([])

    if ylabel is True:
        ax.set_ylabel("Power [dB]")
        ax.set_yticks([-40, -30, -20, -10, 0])
        dax.set_ylabel(r"$\Delta$P [dB]")
    else:
        dax.set_yticklabels([])
        ax.set_yticklabels([])

    # dax = fig.add_subplot(2,1,2)
    dax.scatter(zen_angle, delta_pow, marker=".", s=36, alpha=0.9, color="#63b7af")
    dax.plot(zen_angle, pow_fit, linewidth=2.1, alpha=0.9, color="#ff8264")
    dax.set_ylim([-10, 10])
    # dax.set_xticklabels([])
    dax.set_xlim([-92, 92])

    return ax


def beam_maps(f):
    """Returns pairs of [NS,EW] slices of maps, for each pointing"""
    f_name, _ = Path(f).name.split(".")
    t_name, r_name, _, _ = f_name.split("_")

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

        # rotate maps so slices can be taken
        fee_r = rotate(nside, angle=-np.pi / 4, healpix_array=fee)
        tile_r = rotate(nside, angle=-np.pi / 4, healpix_array=tile)

        # slice the tile and fee maps along NS, EW
        # zenith angle thresh of 70 to determine fit gain factor
        NS_f, EW_f = slice_map(fee_r, 70)
        NS_t, EW_t = tile_map_slice(tile_r, 70)

        gain_NS = fit_gain(map_data=NS_t[0], map_error=NS_t[1], beam=NS_f[0])
        gain_EW = fit_gain(map_data=EW_t[0], map_error=EW_t[1], beam=EW_f[0])

        # slice the tile and fee maps along NS, EW.
        # the above gain factor is applied to full beam slices
        NS_fee, EW_fee = slice_map(fee_r, 90)
        NS_tile, EW_tile = tile_map_slice(tile_r, 90)

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


def plt_grid(tile_name, ref_name):
    # This is an Awesome plot

    f_xx = f"{map_dir}/{tile_name}XX_{ref_name}XX_tile_maps.npz"
    f_yy = f"{map_dir}/{tile_name}YY_{ref_name}YY_tile_maps.npz"

    maps_xx = beam_maps(f_xx)
    maps_yy = beam_maps(f_yy)

    nice_fonts = {
        # Use LaTeX to write all text
        # "text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 14,
        "axes.titlesize": 12,
        "font.size": 14,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 9,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "ytick.color": "w",
        "xtick.color": "w",
        "axes.labelcolor": "w",
        "axes.edgecolor": "w",
    }

    plt.rcParams.update(nice_fonts)

    fig1 = plt.figure(figsize=(11, 6))

    # ax1  = fig1.add_axes([0.01, 0.75, 0.29, 0.22])

    ax1 = plt.subplot(2, 3, 1)
    plt_slice(
        fig=fig1,
        sub=[2, 3, 1],
        zen_angle=maps_xx[0][0][0][2],
        map_slice=maps_xx[0][0][2],
        map_error=maps_xx[0][0][0][1],
        model_slice=maps_xx[0][0][1][0],
        delta_pow=maps_xx[0][0][3],
        pow_fit=maps_xx[0][0][4],
        slice_label="Tile NS",
        model_label="FEE NS",
        ylabel=True,
        title=fr"($i$) NS slice of {tile_name}XX [zenith]",
    )

    ax2 = plt.subplot(2, 3, 2)
    plt_slice(
        fig=fig1,
        sub=[2, 3, 2],
        zen_angle=maps_xx[1][0][0][2],
        map_slice=maps_xx[1][0][2],
        map_error=maps_xx[1][0][0][1],
        model_slice=maps_xx[1][0][1][0],
        delta_pow=maps_xx[1][0][3],
        pow_fit=maps_xx[1][0][4],
        slice_label="Tile NS",
        model_label="FEE NS",
        title=fr"($ii$) NS slice of {tile_name}XX [2]",
    )

    ax3 = plt.subplot(2, 3, 3)
    plt_slice(
        fig=fig1,
        sub=[2, 3, 3],
        zen_angle=maps_xx[2][0][0][2],
        map_slice=maps_xx[2][0][2],
        map_error=maps_xx[2][0][0][1],
        model_slice=maps_xx[2][0][1][0],
        delta_pow=maps_xx[2][0][3],
        pow_fit=maps_xx[2][0][4],
        slice_label="Tile NS",
        model_label="FEE NS",
        title=fr"($iii$) NS slice of {tile_name}XX [4]",
    )

    ax4 = plt.subplot(2, 3, 4)
    plt_slice(
        fig=fig1,
        sub=[2, 3, 4],
        zen_angle=maps_xx[0][1][0][2],
        map_slice=maps_xx[0][1][2],
        map_error=maps_xx[0][1][0][1],
        model_slice=maps_xx[0][1][1][0],
        delta_pow=maps_xx[0][1][3],
        pow_fit=maps_xx[0][1][4],
        slice_label="Tile EW",
        model_label="FEE EW",
        ylabel=True,
        xlabel=True,
        title=fr"($iv$) EW slice of {tile_name}XX [zenith]",
    )

    ax5 = plt.subplot(2, 3, 5)
    plt_slice(
        fig=fig1,
        sub=[2, 3, 5],
        zen_angle=maps_xx[1][1][0][2],
        map_slice=maps_xx[1][1][2],
        map_error=maps_xx[1][1][0][1],
        model_slice=maps_xx[1][1][1][0],
        delta_pow=maps_xx[1][1][3],
        pow_fit=maps_xx[1][1][4],
        slice_label="Tile EW",
        model_label="FEE EW",
        xlabel=True,
        title=fr"($v$) EW slice of {tile_name}XX [2]",
    )

    ax6 = plt.subplot(2, 3, 6)
    plt_slice(
        fig=fig1,
        sub=[2, 3, 6],
        zen_angle=maps_xx[2][1][0][2],
        map_slice=maps_xx[2][1][2],
        map_error=maps_xx[2][1][0][1],
        model_slice=maps_xx[2][1][1][0],
        delta_pow=maps_xx[2][1][3],
        pow_fit=maps_xx[2][1][4],
        slice_label="Tile EW",
        model_label="FEE EW",
        xlabel=True,
        title=fr"($vi$) EW slice of {tile_name}XX [4]",
    )

    plt.tight_layout()
    plt.savefig(
        f"{tile_name}_{ref_name}_map_slices.png",
        transparent=True,
        bbox_inches="tight",
        dpi=420,
    )
    plt.close()


if __name__ == "__main__":

    import argparse
    import numpy as np
    from pathlib import Path
    import concurrent.futures
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gs
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.stats import median_absolute_deviation as mad

    import sys

    sys.path.append("../decode_rf_data")
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
        default="../../../outputs/paper_plots/tile_grids/",
        help="Output directory. Default=../../../outputs/paper_plots/tile_grids",
    )
    parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="../../../outputs/tile_maps/tile_maps_norm/",
        help="Tile map directory. Default=../../../outputs/tile_maps/tile_maps_norm/",
    )
    parser.add_argument(
        "--fee_map",
        metavar="\b",
        default="../../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz",
        help="Healpix FEE map of mwa tile. default=../../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz",
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

    # for r in refs:
    #    for t in tiles:
    #        plt_grid(t, r)

    plt_grid("S07", "rf1")
