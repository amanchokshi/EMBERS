from pathlib import Path

import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import jade
from matplotlib import pyplot as plt
from healpy.rotator import euler_matrix_new as euler


jade, _ = jade()

fee_map = "../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz"
nside = 32


def reproject_map(nside, pixel=None, healpix_array=None):
    """Reproject beam map to be centered around healpix pixel instead of zenith.

    :param nside: Healpix nside
    :param pixel: Healpix index of new phase centre
    :param healpix_array: Input healpix array to be rotated

    :returns:
        - Reprojected healpix map

    """

    # coordinated of new phase centre
    theta, phi = hp.pix2ang(nside, pixel)

    # theta phi values of each pixel
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    new_theta = []
    new_phi = []

    for i in θ - theta:
        if i >= 0:
            new_theta.append(i)
        else:
            new_theta.append(-i)

    for j in ɸ - phi:
        if j <= 2 * np.pi:
            new_phi.append(j)
        else:
            new_phi.append(j - 2 * np.pi)

    new_hp_inds = hp.ang2pix(nside, new_theta, ɸ)
    #  new_hp_inds = hp.ang2pix(nside, θ - theta, ɸ - phi)

    return healpix_array[new_hp_inds]


def plot_healpix(
    data_map=None,
    rot=(0, 90, 180),
    fig=None,
    sub=None,
    title=None,
    vmin=-50,
    vmax=0,
    cmap=None,
    cbar=True,
):
    """Yeesh do some healpix magic to plot the thing

    :param data_map: Healpix input map to plot
    :param fig: Figure number to use
    :param sub: Matplotlib subplot syntax
    :param title: Plot title
    :param vmin: Colormap minimum
    :param vmax: Colormap maximum
    :param cmap: Matplotlib :class:`~matplotlib.colors.ListedColormap`
    :param cbar: If True, plot a colorbar

    :returns:
        - Plot of healpix map

    """

    # Disable cryptic healpy warnings. Can't figure out where they originate
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    hp.delgraticules()
    hp.orthview(
        map=data_map,
        coord="E",
        fig=fig,
        half_sky=True,
        rot=rot,
        xsize=1200,
        title=title,
        sub=sub,
        min=vmin,
        max=vmax,
        cmap=cmap,
        notext=True,
        hold=True,
        cbar=cbar,
        return_projected_map=True,
    )

    hp.graticule(dpar=10, coord="E", color="k", alpha=0.7, dmer=45, lw=0.4, ls=":")

    # Altitude grid
    hp.projtext(
        00.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "0",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )
    hp.projtext(
        30.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "30",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )
    hp.projtext(
        60.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "60",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )

    # NSEW
    hp.projtext(
        80.0 * (np.pi / 180.0),
        000.0 * (np.pi / 180.0),
        r"$N  $",
        coord="E",
        color="w",
        fontweight="light",
        verticalalignment="top",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        090.0 * (np.pi / 180.0),
        r"$E  $",
        coord="E",
        color="w",
        fontweight="light",
        horizontalalignment="right",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        180.0 * (np.pi / 180.0),
        r"$S  $",
        coord="E",
        color="w",
        fontweight="light",
        verticalalignment="bottom",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        270.0 * (np.pi / 180.0),
        r"$W  $",
        coord="E",
        color="w",
        fontweight="light",
        horizontalalignment="left",
    )


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


maps = beam_maps(
    "../embers_out/tile_maps/tile_maps/tile_maps_clean/S06XX_rf0XX_tile_maps.npz"
)[0][0]

#  map_rot = reproject_map(nside, 0, healpix_array=maps)

#  rot_map = plot_healpix(data_map=maps, rot=(90, 75, 90))
#  plt.savefig("rot.png")



vec = hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))
eu_mat = euler(0, 10, 0, deg=True)
rot_map = hp.rotator.rotateVector(eu_mat, vec)
new_inds = hp.vec2pix(nside, rot_map[0], rot_map[1], rot_map[2])
plot_healpix(data_map=maps[new_inds])
plt.savefig("rot.png")
