"""
Compare Beams
-------------

compare measured MWA beam slices to FEE beam models

"""

import argparse
from pathlib import Path

from embers.tile_maps.compare_beams import batch_compare_beam


def main():
    """
    Compare measured MWA beam maps with MWA FEE models using the :func:`~embers.tile_maps.compare_beams.batch_compare_beam` function.

    .. code-block:: console

        $ compare_beams --help

    """

    _parser = argparse.ArgumentParser(
        description="""
        Compare measured MWA beams with FEE models
        """
    )

    _parser.add_argument(
        "--nside", metavar="\b", default=32, type=int, help="Healpix nside. Default=32"
    )

    _parser.add_argument(
        "--fee_map",
        metavar="\b",
        default="embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
        help="MWA FEE healpix model created by embers.mwa_utils.mwa_fee. Default: embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
    )

    _parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="embers_out/tile_maps/tile_maps/tile_maps_clean",
        help="Directory with tile_maps_clean, created by embers.tile_maps.tile_maps.mwa_clean_maps. Default: ../../embers_out/tile_maps/tile_maps/tile_maps_clean",
    )

    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/tile_maps/compare_beams",
        help="Dir where beam comparison plots will be saved. Default=./embers_out/tile_maps/compare_beams",
    )

    _parser.add_argument(
        "--max_cores",
        metavar="\b",
        type=int,
        help="Maximum number of cores to be used by this script. By default all core available cores are used",
    )

    _args = _parser.parse_args()
    _nside = _args.nside
    _fee_map = _args.fee_map
    _map_dir = Path(_args.map_dir)
    _out_dir = Path(_args.out_dir)
    _max_cores = _args.max_cores

    batch_compare_beam(_nside, _fee_map, _map_dir, _out_dir, max_cores=_max_cores)

    print(f"Beam comparison plots saved to: {_out_dir}")
