from pathlib import Path

import healpy as hp
import numpy as np
from healpy.rotator import euler_matrix_new as euler
from scipy import interpolate as interp


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


def beam_rbf_interp(beam_map, nside, dest_nside):

    # Input data
    i_pix = np.arange(hp.nside2npix(nside))
    theta, phi = hp.pix2ang(nside, i_pix)
    power = beam_map[:6144]
    mask = np.isnan(power)
    theta = theta[:6144][~mask]
    phi = phi[:6144][~mask]
    power = power[~mask]
    rbfi = interp.Rbf(theta, phi, power, function="cubic")

    t, p = hp.pix2ang(dest_nside, np.arange(int(hp.nside2npix(dest_nside))))

    beam_interp = rbfi(t, p)
    return beam_interp


if __name__ == "__main__":

    fee_map = "../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz"
    map_dir = "../embers_out/tile_maps/tile_maps/tile_maps_clean"
    out_dir = "./interp_maps"

    nside = 32

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

    tile_pairs = []
    for r in refs:
        for t in tiles:
            tile_pairs.append([t, r])

    for pair in tile_pairs:
        for pol in ["XX", "YY"]:
            tile = f"{pair[0]}{pol}_{pair[1]}{pol}"
            f_tile = f"{map_dir}/{tile}_tile_maps.npz"

            pointings = ["0", "2", "4"]

            # load data from map .npz file
            tile_map = np.load(f_tile, allow_pickle=True)
            fee_m = np.load(fee_map, allow_pickle=True)

            for p in pointings:

                print(f"Crunching tile [{tile}] at pointing [{p}]")

                # Make output directory
                Path(f"{out_dir}/{p}").mkdir(parents=True, exist_ok=True)

                tile = tile_map[p]

                if pol == "XX":
                    fee = fee_m[p][0]
                else:
                    fee = fee_m[p][1]

                # healpix meadian map
                tile_med = np.asarray(
                    [(np.nanmedian(j) if j != [] else np.nan) for j in tile]
                )

                beam_interp = beam_rbf_interp(tile_med, 32, 256)

                np.savez_compressed(
                    f"{out_dir}/{p}/{tile}_{p}_N256", healpix=beam_interp
                )
