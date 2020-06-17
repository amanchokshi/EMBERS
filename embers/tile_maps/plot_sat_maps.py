import sys
import numpy as np
import healpy as hp
from scipy.stats import median_absolute_deviation as mad

sys.path.append("../sat_ephemeris")
from sat_ids import norad_ids

from plot_tile_maps import plot_healpix


def sat_maps(sat):

    f = Path(f"{map_dir}/S07XX_rf0XX_sat_maps.npz")
    f_name, _ = f.name.split(".")
    tile, ref, _, _ = f_name.split("_")

    pointings = ["0", "2", "4", "41"]

    # list of all possible satellites
    sat_ids = list(norad_ids.values())

    # load data from map .npz file
    tile_data = np.load(f, allow_pickle=True)
    tile_data = {key: tile_data[key].item() for key in tile_data}
    tile_map = tile_data["mwa_map"]

    for p in pointings:

        fig = plt.figure(figsize=(8, 10))
        fig.suptitle(f"Satellite [{sat}]: {tile}/{ref} @ {p}", fontsize=16)
        tile_sat_med = [
            (np.median(i) if i != [] else np.nan) for i in tile_data["mwa_map"][p][sat]
        ]
        # tile_sat_scaled = np.asarray([(i - np.nanmax(tile_sat_med[:5000])) for i in tile_sat_med])
        plot_healpix(data_map=np.asarray(tile_sat_med), sub=(1, 1, 1), cmap=jade)
        plt.savefig(
            f"{out_dir}/sat_maps/{p}/{sat}_{p}_{tile}_{ref}_passes.png",
            bbox_inches="tight",
        )
        plt.close()


if __name__ == "__main__":

    import argparse
    import matplotlib

    matplotlib.use("Agg")
    from pathlib import Path
    import concurrent.futures
    import matplotlib.pyplot as plt

    import sys

    sys.path.append("../decode_rf_data")
    from colormap import spectral, jade, kelp

    # Custom spectral colormap
    jade, _ = jade()

    parser = argparse.ArgumentParser(
        description="""
        Plot healpix map of reference data
        """
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./../../outputs/tile_maps/",
        help="Output directory. Default=./../../outputs/tile_maps/",
    )
    parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="./../../outputs/tile_maps/tile_maps_raw",
        help="Output directory. Default=./../../outputs/tile_maps/tile_maps_raw",
    )
    parser.add_argument(
        "--nside",
        metavar="\b",
        type=int,
        default=32,
        help="Healpix Nside. Default = 32",
    )

    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    map_dir = Path(args.map_dir)
    nside = args.nside

    # list of all possible satellites
    sat_ids = list(norad_ids.values())

    pointings = ["0", "2", "4", "41"]

    map_files = [item for item in map_dir.glob("*.npz")]

    # Create output directory tree
    for p in pointings:
        Path(f"{out_dir}/sat_maps/{p}/").mkdir(parents=True, exist_ok=True)

    # plot maps for all sats, for only one tile_ref pair
    # S08XX has good time coverage
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(sat_maps, sat_ids)
    # sat_maps(sat_ids[0])
