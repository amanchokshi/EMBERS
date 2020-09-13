"""
TILE MAPS.
----------
"""

import argparse
import concurrent.futures
import itertools
from itertools import repeat
from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
from embers.rf_tools.colormaps import jade
from embers.tile_maps.beam_utils import plot_healpix
from matplotlib import pyplot as plt
from scipy.interpolate import RectSphereBivariateSpline

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


# https://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent


def polyfit2d(x, y, z, order=1):
    ncols = (order + 1) ** 2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order + 1), range(order + 1))
    for k, (i, j) in enumerate(ij):
        G[:, k] = x ** i * y ** j
    m, _, _, _ = np.linalg.lstsq(G, z, rcond=None)
    return m


def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order + 1), range(order + 1))
    z = np.zeros_like(x)
    for a, (i, j) in zip(m, ij):
        z += a * x ** i * y ** j
    return z


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

    # Fit a 3rd order, 2d polynomial
    m = polyfit2d(x, y, z)

    # Evaluate it on a grid...
    nx, ny = 128, 128
    xx, yy = np.meshgrid(
        np.linspace(x.min(), x.max(), nx), np.linspace(y.min(), y.max(), ny)
    )
    zz = polyval2d(xx, yy, m)

    fig, ax = plt.subplots()
    im = ax.imshow(zz, extent=(-1, 1, 1, -1))
    sc = ax.scatter(x, y, c=z, s=7)
    fig.colorbar(im, ax=ax)
    plt.savefig("surf.png")


tile_gradint(f"{map_dir}/S06XX_rf0XX_tile_maps.npz", 60)
