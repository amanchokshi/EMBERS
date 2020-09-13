"""
TILE TILT.
----------
"""

import argparse
import concurrent.futures
from itertools import repeat
from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
from embers.rf_tools.colormaps import jade
from matplotlib import pyplot as plt
from scipy.linalg import lstsq

matplotlib.use("Agg")
jade, _ = jade()

parser = argparse.ArgumentParser(
    description="""
    Tile Maps paper plot
    """
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots/tile_maps/",
    help="Output directory. Default=../embers_out/paper_plots/tile_maps",
)
parser.add_argument(
    "--map_dir",
    metavar="\b",
    default="../embers_out/tile_maps/tile_maps/tile_maps_clean/",
    help="Tile map directory. Default=../embers_out/tile_maps/tile_maps/tile_maps_clean",
)
parser.add_argument(
    "--fee_map",
    metavar="\b",
    default="../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
    help="Healpix FEE map of mwa tile. default=../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
)
parser.add_argument(
    "--nside", metavar="\b", type=int, default=32, help="Healpix Nside. Default = 32",
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

tile_pairs = []
for r in refs:
    for t in tiles:
        tile_pairs.append([t, r])


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


# https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points

def tile_gradint(f, za):

    maps = beam_maps(f)
    res = maps[0][1]

    pixels = np.arange(hp.nside2npix(nside))
    x, y, _ = hp.pix2vec(nside, pixels)

    # Select central section of beam residual to fit
    pb = int(hp.ang2pix(32, np.radians(za), 0))
    res = res[:pb]
    z = res[~np.isnan(res)]
    x = x[:pb][~np.isnan(res)]
    y = y[:pb][~np.isnan(res)]

    tmp_A = []
    tmp_b = []
    for i in range(len(x)):
        tmp_A.append([x[i], y[i], 1])
        tmp_b.append(z[i])
    b = np.matrix(tmp_b).T
    A = np.matrix(tmp_A)
    fit, residual, rnk, s = lstsq(A, b)

    X, Y = np.meshgrid(np.linspace(-1, 1, 128), np.linspace(-1, 1, 128))
    Z = np.zeros(X.shape)
    for r in range(X.shape[0]):
        for c in range(X.shape[1]):
            Z[r, c] = fit[0] * X[r, c] + fit[1] * Y[r, c] + fit[2]
    fig, ax = plt.subplots()
    im = ax.imshow(Z, extent=(-1, 1, 1, -1))
    sc = ax.scatter(x, y, c=z, s=7)
    fig.colorbar(im, ax=ax)

    plt.savefig("surf.png")


tile_gradint(f"{map_dir}/S06XX_rf0XX_tile_maps.npz", 90)
