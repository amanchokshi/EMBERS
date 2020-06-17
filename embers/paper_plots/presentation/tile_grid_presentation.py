import sys
import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly
from mwa_pb.mwa_sweet_spots import all_grid_points

sys.path.append("../../decode_rf_data")
from colormap import spectral, jade, kelp

# Custom spectral colormap
jade, _ = jade()


def plot_healpix(
    data_map=None,
    fig=None,
    sub=None,
    title=None,
    vmin=None,
    vmax=None,
    cmap=None,
    cbar=True,
):
    """Yeesh do some healpix magic to plot the thing"""

    # Disable cryptic healpy warnings. Can't figure out where they originate
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    hp.delgraticules()
    hp.orthview(
        map=data_map,
        coord="E",
        fig=fig,
        half_sky=True,
        rot=(0, 90, 180),
        xsize=1200,
        title=title,
        sub=sub,
        min=vmin,
        max=vmax,
        cmap=cmap,
        notext=True,
        hold=True,
        cbar=cbar,
        return_projected_map=False,
    )

    hp.graticule(dpar=10, coord="E", color="k", alpha=0.7, dmer=45, lw=0.4, ls=":")

    # Altitude grid
    hp.projtext(
        00.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "0",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )
    hp.projtext(
        30.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "30",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )
    hp.projtext(
        60.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "60",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )

    # NSEW
    hp.projtext(
        80.0 * (np.pi / 180.0),
        000.0 * (np.pi / 180.0),
        r"$N  $",
        coord="E",
        color="w",
        fontweight="light",
        verticalalignment="top",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        090.0 * (np.pi / 180.0),
        r"$E  $",
        coord="E",
        color="w",
        fontweight="light",
        horizontalalignment="right",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        180.0 * (np.pi / 180.0),
        r"$S  $",
        coord="E",
        color="w",
        fontweight="light",
        verticalalignment="bottom",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        270.0 * (np.pi / 180.0),
        r"$W  $",
        coord="E",
        color="w",
        fontweight="light",
        horizontalalignment="left",
    )


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

    fig1 = plt.figure(figsize=(11, 6))

    ax1 = fig1.add_axes([0.01, 0.5, 0.29, 0.44])
    plot_healpix(
        data_map=maps_xx[0][0],
        fig=fig1,
        title=fr"($i$) {tile_name}XX [Zenith]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )
    ax2 = fig1.add_axes([0.30, 0.5, 0.29, 0.44])
    plot_healpix(
        data_map=maps_xx[1][0],
        fig=fig1,
        title=fr"($ii$) {tile_name}XX [2]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )
    ax3 = fig1.add_axes([0.60, 0.5, 0.29, 0.44])
    plot_healpix(
        data_map=maps_xx[2][0],
        fig=fig1,
        title=fr"($iii$) {tile_name}XX [4]",
        cmap=jade,
        vmin=-50,
        vmax=0,
        cbar=False,
    )
    ax3 = plt.gca()
    image = ax3.get_images()[0]
    cax1 = fig1.add_axes([0.91, 0.5, 0.015, 0.44])
    cbar1 = fig1.colorbar(image, cax=cax1, label="Power [dB]")

    ax4 = fig1.add_axes([0.01, 0.0, 0.29, 0.44])
    plot_healpix(
        data_map=maps_xx[0][1],
        fig=fig1,
        title=fr"($iv$) {tile_name}XX - FEE [Zenith]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )
    ax5 = fig1.add_axes([0.30, 0.0, 0.29, 0.44])
    plot_healpix(
        data_map=maps_xx[1][1],
        fig=fig1,
        title=fr"($v$) {tile_name}XX - FEE [2]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )
    ax6 = fig1.add_axes([0.60, 0.0, 0.29, 0.44])
    plot_healpix(
        data_map=maps_xx[2][1],
        fig=fig1,
        title=fr"($vi$) {tile_name}XX - FEE [4]",
        cmap="RdYlGn",
        vmin=-5,
        vmax=5,
        cbar=False,
    )
    ax6 = plt.gca()
    image = ax6.get_images()[0]
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

    import sys

    sys.path.append("../decode_rf_data")
    from colormap import spectral

    # Custom spectral colormap
    cmap = spectral()

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
