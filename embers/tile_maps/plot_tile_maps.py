import sys
import numpy as np
import healpy as hp
from mwa_pb.mwa_sweet_spots import all_grid_points
from scipy.stats import median_absolute_deviation as mad


def plot_healpix(data_map=None, sub=None, title=None, vmin=None, vmax=None, cmap=None):
    """Yeesh do some healpix magic to plot the thing"""

    # Healpix plotting script adapted from Dr. Jack Line's code
    # https://github.com/JLBLine/MWA_ORBCOMM

    # Disable cryptic healpy warnings. Can't figure out where they originate
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    if vmin == None:
        if cmap == None:
            half_sky = hp.orthview(
                map=data_map,
                coord="E",
                half_sky=True,
                xsize=400,
                title=title,
                rot=(0, 90, 180),
                sub=sub,
                notext=True,
                return_projected_map=True,
            )
        else:
            half_sky = hp.orthview(
                map=data_map,
                coord="E",
                half_sky=True,
                xsize=400,
                title=title,
                rot=(0, 90, 180),
                sub=sub,
                cmap=cmap,
                notext=True,
                return_projected_map=True,
            )
    else:
        if cmap == None:
            half_sky = hp.orthview(
                map=data_map,
                coord="E",
                half_sky=True,
                xsize=400,
                rot=(0, 90, 180),
                title=title,
                sub=sub,
                min=vmin,
                max=vmax,
                notext=True,
                return_projected_map=True,
            )
        else:
            half_sky = hp.orthview(
                map=data_map,
                coord="E",
                half_sky=True,
                xsize=400,
                rot=(0, 90, 180),
                title=title,
                sub=sub,
                min=vmin,
                max=vmax,
                cmap=cmap,
                notext=True,
                return_projected_map=True,
            )

    hp.graticule(dpar=10, coord="E", color="k", alpha=0.3, dmer=45)

    # Altitude grid
    hp.projtext(00.0 * (np.pi / 180.0), 225.0 * (np.pi / 180), "0", coord="E")
    hp.projtext(30.0 * (np.pi / 180.0), 225.0 * (np.pi / 180), "30", coord="E")
    hp.projtext(60.0 * (np.pi / 180.0), 225.0 * (np.pi / 180), "60", coord="E")

    # Azimuth grid
    hp.projtext(
        75.0 * (np.pi / 180.0),
        000.0 * (np.pi / 180.0),
        r"$0^\circ$",
        coord="E",
        color="w",
        fontsize=10,
        verticalalignment="top",
    )
    hp.projtext(
        75.0 * (np.pi / 180.0),
        090.0 * (np.pi / 180.0),
        r"$90^\circ$",
        coord="E",
        color="w",
        fontsize=10,
        horizontalalignment="right",
    )
    hp.projtext(
        75.0 * (np.pi / 180.0),
        180.0 * (np.pi / 180.0),
        r"$180^\circ$",
        coord="E",
        color="w",
        fontsize=10,
        verticalalignment="bottom",
    )
    hp.projtext(
        75.0 * (np.pi / 180.0),
        270.0 * (np.pi / 180.0),
        r"$270^\circ$",
        coord="E",
        color="w",
        fontsize=10,
        horizontalalignment="left",
    )

    # NSEW
    hp.projtext(
        90.0 * (np.pi / 180.0),
        000.0 * (np.pi / 180.0),
        r"$N  $",
        coord="E",
        color="k",
        fontsize=14,
        verticalalignment="bottom",
    )
    hp.projtext(
        90.0 * (np.pi / 180.0),
        090.0 * (np.pi / 180.0),
        r"$E  $",
        coord="E",
        color="k",
        fontsize=14,
        horizontalalignment="left",
    )
    hp.projtext(
        90.0 * (np.pi / 180.0),
        180.0 * (np.pi / 180.0),
        r"$S  $",
        coord="E",
        color="k",
        fontsize=14,
        verticalalignment="top",
    )
    hp.projtext(
        90.0 * (np.pi / 180.0),
        270.0 * (np.pi / 180.0),
        r"$W  $",
        coord="E",
        color="k",
        fontsize=14,
        horizontalalignment="right",
    )


def plt_good_maps(f):

    f_name, _ = f.name.split(".")
    tile, ref, _, _ = f_name.split("_")

    pointings = ["0", "2", "4", "41"]

    # load data from map .npz file
    tile_data = np.load(f, allow_pickle=True)

    for p in pointings:

        ## find the pointing center in radians
        # pointing_center_az = np.radians(all_grid_points[int(p)][1])
        # pointing_center_za = np.radians(all_grid_points[int(p)][3])
        #
        ## convert it to a healpix vector
        # pointing_vec = hp.ang2vec(pointing_center_za, pointing_center_az)

        ## find all healpix indices within 10 degrees of pointing center
        # ipix_disc = hp.query_disc(nside=nside, vec=pointing_vec, radius=np.radians(10))

        # healpix meadian map
        tile_map_med = np.asarray(
            [(np.nanmedian(i) if i != [] else np.nan) for i in tile_data[p]]
        )

        ## find the max value within 10 degrees of pointing center
        # ipix_max = np.nanmax(tile_map_med[ipix_disc])
        #
        ## scale map such that the above max is set to 0dB
        # tile_map_scaled = np.asarray([(i - ipix_max) for i in tile_map_med])

        fig = plt.figure(figsize=(8, 10))
        fig.suptitle(f"Good Map: {tile}/{ref} @ {p}", fontsize=16)
        plot_healpix(data_map=tile_map_med, sub=(1, 1, 1), cmap=jade, vmin=-50, vmax=0)
        plt.savefig(
            f"{out_dir}/{p}/tile_maps/{tile}_{ref}_{p}_good_map.png",
            bbox_inches="tight",
        )
        plt.close()

        # Plot MAD
        tile_map_mad = []
        for j in tile_data[p]:
            if j != []:
                j = np.asarray(j)
                j = j[~np.isnan(j)]
                tile_map_mad.append(mad(j))
            else:
                tile_map_mad.append(np.nan)

        vmin = np.nanmin(tile_map_mad)
        vmax = np.nanmax(tile_map_mad)

        fig = plt.figure(figsize=(8, 10))
        fig.suptitle(f"Good Map MAD: {tile}/{ref} @ {p}", fontsize=16)
        plot_healpix(
            data_map=np.asarray(tile_map_mad),
            sub=(1, 1, 1),
            cmap=jade,
            vmin=vmin,
            vmax=vmax,
        )
        plt.savefig(
            f"{out_dir}/{p}/tile_errors/{tile}_{ref}_{p}_good_map_errors.png",
            bbox_inches="tight",
        )
        plt.close()

        # Plot counts in pix
        tile_map_counts = [len(np.array(i)[~np.isnan(i)]) for i in tile_data[p]]

        fig = plt.figure(figsize=(8, 10))
        fig.suptitle(f"Good Map Counts: {tile}/{ref} @ {p}", fontsize=16)
        plot_healpix(
            data_map=np.asarray(tile_map_counts),
            sub=(1, 1, 1),
            cmap=jade,
            vmin=0,
            vmax=80,
        )
        plt.savefig(
            f"{out_dir}/{p}/tile_counts/{tile}_{ref}_{p}_good_map_counts.png",
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
        Plot healpix maps of normalized tile maps
        """
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./../../outputs/tile_maps/good_maps",
        help="Output directory. Default=./../../outputs/tile_maps/good_maps",
    )
    parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="./../../outputs/tile_maps/tile_maps_norm",
        help="Output directory. Default=./../../outputs/tile_maps/tile_maps_norm",
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

    pointings = ["0", "2", "4", "41"]

    map_files = [item for item in map_dir.glob("*.npz")]

    # Create output directory tree
    for p in pointings:
        Path(f"{out_dir}/{p}/tile_maps/").mkdir(parents=True, exist_ok=True)
        Path(f"{out_dir}/{p}/tile_counts/").mkdir(parents=True, exist_ok=True)
        Path(f"{out_dir}/{p}/tile_errors/").mkdir(parents=True, exist_ok=True)

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(plt_good_maps, map_files)
