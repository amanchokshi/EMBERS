import argparse
import concurrent.futures
from itertools import repeat
from pathlib import Path

import matplotlib
import numpy as np
from embers.rf_tools.colormaps import spectral
from embers.tile_maps.beam_utils import (chisq_fit_gain,
                                         healpix_cardinal_slices, map_slices,
                                         poly_fit, rotate_map)
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

matplotlib.use("Agg")
_spec, _ = spectral()

parser = argparse.ArgumentParser(
    description="""
    Tile Slices paper plot
    """
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots/tile_slices/",
    help="Output directory. Default=../embers_out/paper_plots/tile_slices",
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

refs = ["rf1"]

tile_pairs = []
for r in refs:
    for t in tiles:
        tile_pairs.append([t, r])


def beam_slices(map_file, fee_map, nside):
    """Returns pairs of [NS,EW] slices of maps, for each pointing"""

    t_name, r_name, _, _ = Path(map_file).stem.split("_")

    pointings = ["0", "2", "4"]

    maps = []

    # load data from map .npz file
    tile_map = np.load(map_file, allow_pickle=True)
    fee_m = np.load(fee_map, allow_pickle=True)

    for p in pointings:

        tile = tile_map[p]

        if "XX" in t_name:
            fee = fee_m[p][0]
        else:
            fee = fee_m[p][1]

        # rotate maps so slices can be taken
        fee_r = rotate_map(nside, angle=-np.pi / 4, healpix_array=fee)
        tile_r = rotate_map(nside, angle=-np.pi / 4, healpix_array=tile)

        # slice the tile and fee maps along NS, EW
        # zenith angle thresh of 70 to determine fit gain factor
        NS_f, EW_f = healpix_cardinal_slices(nside, fee_r, 70)
        NS_t, EW_t = map_slices(nside, tile_r, 70)

        gain_NS = chisq_fit_gain(data=NS_t[0], model=NS_f[0])
        gain_EW = chisq_fit_gain(data=EW_t[0], model=EW_f[0])

        # slice the tile and fee maps along NS, EW.
        # the above gain factor is applied to full beam slices
        NS_fee, EW_fee = healpix_cardinal_slices(nside, fee_r, 90)
        NS_tile, EW_tile = map_slices(nside, tile_r, 90)

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


slice_0_NS = {}
slice_0_EW = {}
slice_2_NS = {}
slice_2_EW = {}
slice_4_NS = {}
slice_4_EW = {}

for tile in tile_pairs:

    tile_name = tile[0]
    ref_name = tile[1]

    for pol in ["XX", "YY"]:
        f_tile = f"{map_dir}/{tile_name}{pol}_{ref_name}{pol}_tile_maps.npz"
        try:
            slices = beam_slices(f_tile, fee_map, nside)
            slice_0_NS[f"{tile_name}{pol}"] = slices[0][0]
            slice_0_EW[f"{tile_name}{pol}"] = slices[0][1]
            slice_2_NS[f"{tile_name}{pol}"] = slices[1][0]
            slice_2_EW[f"{tile_name}{pol}"] = slices[1][1]
            slice_4_NS[f"{tile_name}{pol}"] = slices[2][0]
            slice_4_EW[f"{tile_name}{pol}"] = slices[2][1]
        except Exception as e:
            print(e)


plt.style.use("seaborn")

nice_fonts = {
    "font.family": "sans-serif",
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "font.size": 8,
    "legend.fontsize": 6,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
}


plt.rcParams.update(nice_fonts)

ax = plt.figure(figsize=(6, 5))

tile_keys = list(slice_0_NS.keys())
colors = _spec(np.linspace(0.17, 0.9, len(tile_keys)))

plt.plot(
    slice_4_NS["S06YY"][0][2],
    slice_4_NS["S06YY"][1][0],
    color="black",
    lw=2,
    label="FEE",
)

for i, k in enumerate(tile_keys):
    if "YY" in k:
        plt.scatter(
            slice_4_NS[k][0][2],
            slice_4_NS[k][2],
            label=k,
            s=21,
            color=colors[i],
            alpha=0.88,
            edgecolor="black",
            linewidth=0.2,
        )

leg = ax.legend(frameon=True, handlelength=1)
leg.get_frame().set_facecolor("white")
for le in leg.legendHandles:
    le.set_alpha(1)
plt.xlabel("Zenith Angle [deg]")
plt.ylabel("Power [dB]")
plt.ylim([-55, 3])
plt.tight_layout()
plt.savefig(f"{out_dir}/YY_4.pdf")
