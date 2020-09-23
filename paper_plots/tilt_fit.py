from pathlib import Path

import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import jade
from embers.tile_maps.beam_utils import plot_healpix
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

    t, p = hp.pix2ang(128, np.arange(int(hp.nside2npix(128))))

    p_128 = rbfi(t, p)
    #  plot_healpix(data_map=p_128, vmin=-50, vmax=0)
    #  plt.savefig("rbf.png")
    return p_128


b_128 = beam_rbf(b_map)
f_128 = beam_rbf(fee)

res_128 = b_128 - f_128
mask = np.where(f_128 < -30)
res_128[mask] = np.nan

plot_healpix(data_map=res_128, vmin=-5, vmax=5)
plt.savefig("res_rbf.png")
