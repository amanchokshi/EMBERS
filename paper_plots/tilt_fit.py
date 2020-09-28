# Phase beam maps to be centered around sucessive healpix pixels, and compute mean residuals
# Make pix vs mean residual plot and save data dictionary to json file

import json
from pathlib import Path

import healpy as hp
import numpy as np
from healpy.rotator import euler_matrix_new as euler
from matplotlib import pyplot as plt


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


def pointing_residual(nside, pixel, beam=None, fee=None):

    # Select primary and secondary lobes
    mask = np.where(fee < -30)

    rot_beam = reproject_map(nside, pixel, healpix_array=beam)

    resi = rot_beam - fee
    resi[mask] = np.nan

    # Select central 40 deg of beam. Primary + edges of secondary lobes
    ind = int(hp.ang2pix(nside, np.radians(40), 0))
    resi = resi[:ind]

    # Mean of residuals as proxy for pointing errors
    resi_mean = np.nanmean(resi)

    return resi_mean


if __name__ == "__main__":

    out_dir = "./interp_maps"

    nside = 256

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

    # Make output directory
    Path(f"{out_dir}/tile_tilt").mkdir(parents=True, exist_ok=True)

    pix_residuals = {}

    for pair in tile_pairs:
        for pol in ["XX", "YY"]:

            pointings = ["0", "2", "4"]

            for p in pointings:

                try:

                    tile = f"{pair[0]}{pol}_{pair[1]}{pol}"
                    beam = f"{out_dir}/{p}/{tile}_{p}_N256.npz"
                    fee = f"{out_dir}/FEE/fee_{pol}_{p}_N256.npz"
                    beam = np.load(beam)["healpix"]
                    fee = np.load(fee)["healpix"]

                    pix_res = []
                    for pix in range(1500):
                        res_mean = pointing_residual(nside, pix, beam=beam, fee=fee)
                        pix_res.append(res_mean)

                    plt.style.use("seaborn")
                    plt.scatter(range(1500), pix_res, color="seagreen", s=7)
                    plt.tight_layout()
                    plt.savefig(
                        f"{out_dir}/tile_tilt/{pair[0]}{pol}_{pair[1]}{pol}_{p}res.png"
                    )
                    plt.close()
                    print(f"{pair[0]}{pol}_{pair[1]}{pol}_{p}res.png saved")
                    pix_residuals[f"{pair[0]}{pol}_{pair[1]}{pol}_{p}"] = pix_res

                except Exception as e:
                    print(e)

    with open(f"{out_dir}/pix_resi.json", "w") as f:
        json.dump(pix_residuals, f, indent=4)
