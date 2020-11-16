"""
Null Test
---------

perform null tests on reference rf data and reference beam models

"""

import argparse
from pathlib import Path

from embers.tile_maps.null_test import null_test


def main():
    """
    Perform a null test of the reference antennas using the :func:`~embers.tile_maps.null_test.null_test` function.

    .. code-block:: console

        $ null_test --help

    """

    _parser = argparse.ArgumentParser(
        description="""
        Create MWA FEE beam models at multiple pointings with dipoles flagged.
        """
    )

    _parser.add_argument(
        "--nside", metavar="\b", default=32, type=int, help="Healpix nside. Default=32"
    )

    _parser.add_argument(
        "--za_max",
        metavar="\b",
        default=90,
        type=int,
        help="Maximum zenith angle upto which to perform the null test. Default: 90 deg",
    )

    _parser.add_argument(
        "--ref_model",
        metavar="\b",
        default="embers_out/tile_maps/ref_models/ref_dipole_models.npz",
        help="Reference feko healpix model created by embers.tile_maps.ref_fee_healpix. Default: embers_out/tile_maps/ref_models/ref_dipole_models.npz",
    )

    _parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="embers_out/tile_maps/tile_maps/tile_maps_raw",
        help="Directory with tile_maps_raw, created by embers.tile_maps.tile_maps.project_tile_healpix. Default: embers_out/tile_maps/tile_maps/tile_maps_raw",
    )

    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/tile_maps/null_test",
        help="Dir where null tests will be saved. Default=./embers_out/tile_maps/null_test",
    )

    _args = _parser.parse_args()
    _nside = _args.nside
    _za_max = _args.za_max
    _ref_model = _args.ref_model
    _map_dir = Path(_args.map_dir)
    _out_dir = Path(_args.out_dir)

    print(f"Null tests saved to {_out_dir}")
    null_test(_nside, _za_max, _ref_model, _map_dir, _out_dir)
