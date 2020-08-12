"""
NULL TEST.
----------

"""

import argparse
from pathlib import Path

from embers.tile_maps.null_test import null_test

parser = argparse.ArgumentParser(
    description="""
    Null Test paper plot
    """
)

parser.add_argument(
    "--nside", metavar="\b", default=32, type=int, help="Healpix nside. Default=32"
)

parser.add_argument(
    "--za_max",
    metavar="\b",
    default=80,
    type=int,
    help="Maximum zenith angle upto which to perform the null test. Default: 80 deg",
)

parser.add_argument(
    "--ref_model",
    metavar="\b",
    default="../embers_out/tile_maps/ref_models/ref_dipole_models.npz",
    help="Reference feko healpix model created by embers.tile_maps.ref_fee_healpix. Default: ../embers_out/tile_maps/ref_models/ref_dipole_models.npz",
)

parser.add_argument(
    "--map_dir",
    metavar="\b",
    default="../embers_out/tile_maps/tile_maps/tile_maps_raw",
    help="Directory with tile_maps_raw, created by embers.tile_maps.tile_maps.project_tile_healpix. Default: ../embers_out/tile_maps/tile_maps/tile_maps_raw",
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="../embers_out/paper_plots",
    help="Output directory. Default=../embers_out/paper_plots",
)

args = parser.parse_args()
nside = args.nside
za_max = args.za_max
ref_model = args.ref_model
map_dir = Path(args.map_dir)
out_dir = Path(args.out_dir)


null_test(nside, za_max, ref_model, map_dir, out_dir)
print(f"NULL TEST saved to {out_dir}")
