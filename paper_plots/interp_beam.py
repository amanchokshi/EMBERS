# Interpolate healpix map using a radial basis function
# N=33 -> N=265. Uses a massive amount of memory (All ram + ~60GB of swap)
# Takes 18 hours to complete for all maps

from pathlib import Path

import healpy as hp
import numpy as np
from scipy import interpolate as interp


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

    # Interpolate FEE to nside = 256
    fee_m = np.load(fee_map, allow_pickle=True)
    # Make output directory
    Path(f"{out_dir}/FEE").mkdir(parents=True, exist_ok=True)

    for p in ["0", "2", "4"]:

        for pol in ["XX", "YY"]:

            print(f"Crunching FEE [{pol}] at pointing [{p}]")

            if pol == "XX":
                fee = fee_m[p][0]
            else:
                fee = fee_m[p][1]

            beam_interp = beam_rbf_interp(fee, 32, 256)

            np.savez_compressed(
                f"{out_dir}/FEE/fee_{pol}_{p}_N256", healpix=beam_interp
            )

    # Interpolate beams to nside = 256
    for pair in tile_pairs:
        for pol in ["XX", "YY"]:
            tile = f"{pair[0]}{pol}_{pair[1]}{pol}"
            f_tile = f"{map_dir}/{tile}_tile_maps.npz"

            pointings = ["0", "2", "4"]

            # load data from map .npz file
            tile_map = np.load(f_tile, allow_pickle=True)

            for p in pointings:

                print(f"Crunching tile [{tile}] at pointing [{p}]")

                # Make output directory
                Path(f"{out_dir}/{p}").mkdir(parents=True, exist_ok=True)

                tile_p = tile_map[p]

                # healpix meadian map
                tile_med = np.asarray(
                    [(np.nanmedian(j) if j != [] else np.nan) for j in tile_p]
                )

                beam_interp = beam_rbf_interp(tile_med, 32, 256)

                np.savez_compressed(
                    f"{out_dir}/{p}/{tile}_{p}_N256", healpix=beam_interp
                )
