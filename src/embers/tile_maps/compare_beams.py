"""
Compare Beams
-------------

A set of tools to compare measured MWA beam maps with FEE models

"""

import concurrent.futures
from itertools import repeat
from pathlib import Path

import matplotlib
import numpy as np
from embers.rf_tools.colormaps import jade
from embers.tile_maps.beam_utils import (chisq_fit_gain,
                                         healpix_cardinal_slices, map_slices,
                                         plot_healpix, plt_slice, poly_fit,
                                         rotate_map)
from matplotlib import pyplot as plt

matplotlib.use("Agg")
jade, _ = jade()


def beam_slice(nside, tile_map, fee_map, out_dir):
    """Compare slices of measured beam maps and FEE models.

    NS & EW slices of measured MWA beam maps are compared to corresponding slices of FEE models. Complete sky maps are plotted to display
    power gradients across the beam.

    :param nside: Healpix nside
    :param tile_map: Clean MWA tile map created by :func:`~embers.tile_maps.tile_maps.mwa_clean_maps`
    :param fee_map: MWA FEE model created  my :func:`~embers.mwa_utils.mwa_fee`
    :param out_dir: Path to output directory where diagnostic plots will be saved
    """

    t_name, r_name, _, _ = tile_map.stem.split("_")

    pointings = ["0", "2", "4", "41"]

    # load data from map .npz file
    tile_map = np.load(tile_map, allow_pickle=True)
    fee_m = np.load(fee_map, allow_pickle=True)

    # MWA beam pointings
    pointings = ["0", "2", "4", "41"]

    for p in pointings:

        Path(f"{out_dir}/{p}/").mkdir(parents=True, exist_ok=True)

        try:
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

            # Visualize the tile map and diff map
            # healpix meadian map
            tile_med = np.asarray(
                [(np.nanmedian(j) if j != [] else np.nan) for j in tile]
            )

            residuals = tile_med - fee
            residuals[np.where(fee < -30)] = np.nan
            residuals[np.where(tile_med == np.nan)] = np.nan

            # This is an Awesome plot
            plt.style.use("seaborn")
            fig1 = plt.figure(figsize=(10, 8))
            ax = plt.gca()
            ax.set_axis_off()

            plt.gca()
            plt_slice(
                fig=fig1,
                sub=(2, 2, 1),
                zen_angle=NS_tile[2],
                map_slice=NS_tile_med,
                map_error=NS_tile[1],
                model_slice=NS_fee[0],
                delta_pow=del_NS,
                pow_fit=fit_NS,
                slice_label="Tile NS",
                model_label="FEE NS",
                xlim=[-90, 90],
                ylim=[-54, 4],
            )

            fig1.add_axes([0.48, 0.52, 0.48, 0.43])
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
            ax1 = plt.gca()
            image = ax1.get_images()[0]
            cax = fig1.add_axes([0.92, 0.52, 0.015, 0.43])
            fig1.colorbar(image, cax=cax, label="dB")

            plt.gca()
            plt_slice(
                fig=fig1,
                sub=(2, 2, 3),
                zen_angle=EW_tile[2],
                map_slice=EW_tile_med,
                map_error=EW_tile[1],
                model_slice=EW_fee[0],
                delta_pow=del_EW,
                pow_fit=fit_EW,
                slice_label="Tile EW",
                model_label="FEE EW",
                xlabel=True,
                xlim=[-90, 90],
                ylim=[-54, 4],
            )

            fig1.add_axes([0.48, 0.02, 0.48, 0.43])
            plot_healpix(
                data_map=residuals,
                sub=(2, 2, 4),
                fig=fig1,
                title="diff map",
                cmap="RdYlGn",
                vmin=-10,
                vmax=5,
                cbar=False,
            )
            ax2 = plt.gca()
            image = ax2.get_images()[0]
            cax = fig1.add_axes([0.92, 0.02, 0.015, 0.43])
            fig1.colorbar(image, cax=cax, label="dB")

            plt.tight_layout()
            plt.savefig(f"{out_dir}/{p}/{t_name}_{r_name}_{p}_beam_slices.png")
            plt.close()

        except Exception as e:
            print(e)


def batch_compare_beam(nside, fee_map, map_dir, out_dir, max_cores=None):
    """Batch compare multiple beam maps

    :param nside: Healpix nside
    :param fee_map: MWA FEE model created  my :func:`~embers.mwa_utils.mwa_fee`
    :param map_dir: Path to dir with clean MWA tile maps created by :func:`~embers.tile_maps.tile_maps.mwa_clean_maps`
    :param out_dir: Path to output directory where diagnostic plots will be saved
    :param max_cores: Maximum number of cores to be used by this script. Default=None, which means that all available cores are used
    """

    # make output dir if it doesn't exist
    out_dir.mkdir(parents=True, exist_ok=True)

    # find all map files
    map_files = [item for item in map_dir.glob("*.npz")]

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_cores) as executor:
        executor.map(
            beam_slice, repeat(nside), map_files, repeat(fee_map), repeat(out_dir)
        )
