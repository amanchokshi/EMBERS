import sys
import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly
from mwa_pb.mwa_sweet_spots import all_grid_points

sys.path.append("../decode_rf_data")
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
    sub=None,
    zen_angle=None,
    map_slice=None,
    map_error=None,
    model_slice=None,
    delta_pow=None,
    pow_fit=None,
    slice_label=None,
    model_label=None,
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

    ax.set_ylabel("Power (dB)")
    # ax.set_ylim(bottom=-30)
    ax.set_xlabel("Zenith Angle (degrees)")
    ax.legend()
    ax.set_xlim([-92, 92])
    ax.set_ylim([-50, 5])

    divider = make_axes_locatable(ax)
    dax = divider.append_axes("bottom", size="30%", pad=0.1)

    # dax = fig.add_subplot(2,1,2)
    dax.scatter(zen_angle, delta_pow, marker=".", color="#27296d")
    dax.plot(zen_angle, pow_fit, linewidth=1.2, alpha=0.9, color="#ff8264")
    dax.set_ylabel("Data - Model (dB)")
    dax.set_ylim([-10, 10])
    # dax.set_xticklabels([])
    dax.set_xlim([-92, 92])

    return ax


def plot_healpix(
    data_map=None,
    fig=None,
    sub=None,
    title=None,
    vmin=None,
    vmax=None,
    cmap=None,
    cbar=True,
):
    """Yeesh do some healpix magic to plot the thing"""

    # Disable cryptic healpy warnings. Can't figure out where they originate
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    hp.orthview(
        map=data_map,
        coord="E",
        fig=fig,
        half_sky=True,
        rot=(0, 90, 180),
        xsize=600,
        title=title,
        sub=sub,
        min=vmin,
        max=vmax,
        cmap=cmap,
        notext=True,
        hold=True,
        cbar=cbar,
    )

    hp.graticule(dpar=10, coord="E", color="k", alpha=0.3, dmer=45)

    # Altitude grid
    hp.projtext(
        00.0 * (np.pi / 180.0), 225.0 * (np.pi / 180), "0", color="w", coord="E"
    )
    hp.projtext(
        30.0 * (np.pi / 180.0), 225.0 * (np.pi / 180), "30", color="w", coord="E"
    )
    hp.projtext(
        60.0 * (np.pi / 180.0), 225.0 * (np.pi / 180), "60", color="w", coord="E"
    )

    # NSEW
    hp.projtext(
        80.0 * (np.pi / 180.0),
        000.0 * (np.pi / 180.0),
        r"$N  $",
        coord="E",
        color="w",
        verticalalignment="top",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        090.0 * (np.pi / 180.0),
        r"$E  $",
        coord="E",
        color="w",
        horizontalalignment="right",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        180.0 * (np.pi / 180.0),
        r"$S  $",
        coord="E",
        color="w",
        verticalalignment="bottom",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        270.0 * (np.pi / 180.0),
        r"$W  $",
        coord="E",
        color="w",
        horizontalalignment="left",
    )


def beam_slice(f):

    f_name, _ = f.name.split(".")
    t_name, r_name, _, _ = f_name.split("_")

    pointings = ["0", "2", "4", "41"]

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

        # Visualize the tile map and diff map
        # healpix meadian map
        tile_med = np.asarray([(np.nanmedian(j) if j != [] else np.nan) for j in tile])

        residuals = tile_med - fee
        residuals[np.where(fee < -30)] = np.nan
        residuals[np.where(tile_med == np.nan)] = np.nan

        # This is an Awesome plot
        plt.style.use("seaborn")
        fig1 = plt.figure(figsize=(10, 8))

        ax1 = plt.subplot(2, 2, 1)
        plt_slice(
            fig=fig1,
            sub=221,
            zen_angle=NS_tile[2],
            map_slice=NS_tile_med,
            map_error=NS_tile[1],
            model_slice=NS_fee[0],
            delta_pow=del_NS,
            pow_fit=fit_NS,
            slice_label="Tile NS",
            model_label="FEE NS",
        )

        ax2 = fig1.add_axes([0.48, 0.52, 0.48, 0.43])
        plot_healpix(
            data_map=tile_med,
            sub=(2, 2, 2),
            fig=fig1,
            title="tile map",
            cmap=jade,
            vmin=-50,
            vmax=0,
            cbar=False,
        )
        ax7 = plt.gca()
        image = ax7.get_images()[0]
        cax = fig1.add_axes([0.92, 0.52, 0.015, 0.43])
        cbar = fig1.colorbar(image, cax=cax, label="dB")

        ax3 = plt.subplot(2, 2, 3)
        plt_slice(
            fig=fig1,
            sub=223,
            zen_angle=EW_tile[2],
            map_slice=EW_tile_med,
            map_error=EW_tile[1],
            model_slice=EW_fee[0],
            delta_pow=del_EW,
            pow_fit=fit_EW,
            slice_label="Tile EW",
            model_label="FEE EW",
        )

        ax4 = fig1.add_axes([0.48, 0.02, 0.48, 0.43])
        plot_healpix(
            data_map=residuals,
            sub=(2, 2, 4),
            fig=fig1,
            title="diff map",
            cmap="inferno",
            vmin=-10,
            vmax=5,
            cbar=False,
        )
        ax8 = plt.gca()
        image = ax8.get_images()[0]
        cax = fig1.add_axes([0.92, 0.02, 0.015, 0.43])
        cbar = fig1.colorbar(image, cax=cax, label="dB")

        plt.tight_layout()
        plt.savefig(f"{out_dir}/{p}/{t_name}_{r_name}_{p}_beam_slices.png")
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
        default="../../outputs/tile_maps/compare_beams/",
        help="Output directory. Default=../../outputs/tile_maps/compare_beams/",
    )
    parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="../../outputs/tile_maps/tile_maps_norm/",
        help="Tile map directory. Default=../../outputs/tile_maps/tile_maps_norm/",
    )
    parser.add_argument(
        "--fee_map",
        metavar="\b",
        default="../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz",
        help="Healpix FEE map of mwa tile. default=../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz",
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

    # MWA beam pointings
    pointings = ["0", "2", "4", "41"]

    # Create output directory structure
    for p in pointings:
        Path(f"{out_dir}/{p}/").mkdir(parents=True, exist_ok=True)

    # find all map files
    map_files = [item for item in map_dir.glob("*.npz")]

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(beam_slice, map_files)
#    beam_slice(map_files[0])
