from pathlib import Path

import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import jade
from embers.tile_maps.beam_utils import plot_healpix
from healpy.rotator import euler_matrix_new as euler
from matplotlib import pyplot as plt
from scipy import interpolate as interp

jade, _ = jade()

fee_map = "../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz"
nside = 32


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


def reproject_map(nside, pixel, healpix_array=None):
    """Reproject beam map to be centered around healpix pixel instead of zenith.
    :param nside: Healpix nside
    :param pixel: Healpix index of new phase centre
    :param healpix_array: Input healpix array to be rotated
    :returns:
        - Reprojected healpix map
    """

    # coordinated of new phase centre
    #  theta, phi = hp.pix2ang(nside, pixel)
    # Theta moves primary beam eastward
    # Phi rotates primary beam anticlockwise. Or, E --> N
    #  theta, phi = [1, -45]
    theta, phi = hp.pix2ang(nside, pixel)

    vec = hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))
    eu_mat = euler(-phi, -theta, phi, deg=False)
    rot_map = hp.rotator.rotateVector(eu_mat, vec)
    new_hp_inds = hp.vec2pix(nside, rot_map[0], rot_map[1], rot_map[2])

    return healpix_array[new_hp_inds]


b_map = beam_maps(
    "../embers_out/tile_maps/tile_maps/tile_maps_clean/S06XX_rf0XX_tile_maps.npz"
)[0][0]

fee_m = np.load(fee_map, allow_pickle=True)
fee = fee_m["0"][0]


def beam_rbf(beam_map):

    # Input data
    i_pix = np.arange(hp.nside2npix(nside))
    theta, phi = hp.pix2ang(nside, i_pix)
    power = beam_map[:6144]
    mask = np.isnan(power)
    theta = theta[:6144][~mask]
    phi = phi[:6144][~mask]
    power = power[~mask]
    rbfi = interp.Rbf(theta, phi, power, function="cubic")

    t, p = hp.pix2ang(512, np.arange(int(hp.nside2npix(512))))

    p_128 = rbfi(t, p)
    #  plot_healpix(data_map=p_128, vmin=-50, vmax=0)
    #  plt.savefig("rbf.png")
    return p_128


b_128 = beam_rbf(b_map)
f_128 = beam_rbf(fee)

for i in range(141):
    b_128 = reproject_map(128, i, b_128)
    res_128 = b_128 - f_128
    mask = np.where(f_128 < -30)
    res_128[mask] = np.nan

    plot_healpix(data_map=res_128, vmin=-5, vmax=5)
    plt.savefig(f"res_rot_{i}.png")
    plt.close()
