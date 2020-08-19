import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly
from mwa_pb.mwa_sweet_spots import all_grid_points

from embers.rf_tools.colormaps import jade
from embers.tile_maps.beam_utils import plot_healpix

# Custom spectral colormap
jade, _ = jade()

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


def plt_grid(tile_name, ref_name):
    # This is an Awesome plot

    f_xx = f"{map_dir}/{tile_name}XX_{ref_name}XX_tile_maps.npz"
    f_yy = f"{map_dir}/{tile_name}YY_{ref_name}YY_tile_maps.npz"

    maps_xx = beam_maps(f_xx)
    maps_yy = beam_maps(f_yy)

    nice_fonts = {
        # Use LaTeX to write all text
        # "text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 14,
        "axes.titlesize": 14,
        "font.size": 14,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10  # ,
        # "ytick.color" : "w",
        # "xtick.color" : "w",
        # "axes.labelcolor" : "w",
        # "axes.edgecolor" : "w"
    }
    plt.rcParams.update(nice_fonts)
    # plt.rcParams['grid.color'] = '#cccccc'
    # plt.rcParams['grid.linestyle'] = ':'
    # plt.rcParams['grid.linewidth'] = 0.3
    # plt.style.use('seaborn')

    fig1 = plt.figure(figsize=(6, 11))

    ax1 = fig1.add_axes([0.01, 0.5, 0.95, 0.44])
    plot_healpix(
        data_map=maps_xx[0][0],
        fig=fig1,
        title=fr"($i$) {tile_name}XX [Zenith]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )
    ax1 = plt.gca()
    image = ax1.get_images()[0]
    cax1 = fig1.add_axes([0.91, 0.5, 0.015, 0.44])
    cbar1 = fig1.colorbar(image, cax=cax1, label="Power [dB]")

    ax2 = fig1.add_axes([0.01, 0.0, 0.95, 0.44])
    plot_healpix(
        data_map=maps_xx[0][1],
        fig=fig1,
        title=fr"($iv$) {tile_name}XX - FEE [Zenith]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )
    ax2 = plt.gca()
    image = ax2.get_images()[0]
    cax2 = fig1.add_axes([0.91, 0.0, 0.015, 0.44])
    cbar2 = fig1.colorbar(image, cax=cax2, label="Power [dB]")

    plt.savefig(f"{tile_name}_{ref_name}_maps.png", bbox_inches="tight", dpi=420)
    plt.close()


if __name__ == "__main__":

    import argparse
    import numpy as np
    from pathlib import Path
    import concurrent.futures
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gs
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.stats import median_absolute_deviation as mad


    parser = argparse.ArgumentParser(
        description="""
        Plot healpix map of reference data
        """
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="../../../outputs/paper_plots/tile_grids/",
        help="Output directory. Default=../../../outputs/paper_plots/tile_grids",
    )
    parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="../../../outputs/tile_maps/tile_maps_norm/",
        help="Tile map directory. Default=../../../outputs/tile_maps/tile_maps_norm/",
    )
    parser.add_argument(
        "--fee_map",
        metavar="\b",
        default="../../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz",
        help="Healpix FEE map of mwa tile. default=../../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz",
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

    # for r in refs:
    #    for t in tiles:
    #        plt_grid(t, r)

    plt_grid("S07", "rf1")
