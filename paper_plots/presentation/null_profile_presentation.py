import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly


def ring_indices(nside=None):
    """Healpix pix indices of NS, EW slices and above horizon"""
    # theta phi values of each pixel
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    za = np.arange(0, 87, 2)

    ring_indices = []
    for i in range(len(za)):
        if i < 43:
            r_indices = np.where(
                np.logical_and(θ >= np.radians(za[i]), θ < np.radians(za[i + 1]))
            )[0]
            ring_indices.append(r_indices)
    za_cen = za[:-1] + 1
    return [za_cen, ring_indices]


def good_maps(ref_map):
    """Creates a ref map with only good satellites"""

    pointings = ["0", "2", "4"]

    # load data from map .npz file
    f = Path(f"{map_dir}/{ref_map}")
    tile_data = np.load(f, allow_pickle=True)
    tile_data = {key: tile_data[key].item() for key in tile_data}
    ref_map = tile_data["ref_map"]

    # Good sats from which to make plots
    good_sats = [
        25338,
        25984,
        25985,
        28654,
        40086,
        40087,
        40091,
        41179,
        41180,
        41182,
        41183,
        41184,
        41185,
        41187,
        41188,
        41189,
        44387,
    ]

    # Empty good map
    good_map = [[] for pixel in range(hp.nside2npix(nside))]

    for p in pointings:

        # append to good map from all good sat data
        for sat in good_sats:
            for pix in range(hp.nside2npix(nside)):
                good_map[pix].extend(ref_map[p][sat][pix])

    good_map = [np.nanmedian(pixel) for pixel in good_map]

    return good_map


def poly_fit(x, y, order):
    """Fit polynominal of order to data"""

    x = np.asarray(x)
    y = np.asarray(y)

    coefs = poly.polyfit(x, y, order)
    fit = poly.polyval(x, coefs)
    return fit


if __name__ == "__main__":

    import argparse
    import numpy as np
    from pathlib import Path
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gs
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.stats import median_absolute_deviation as mad

    import sys

    sys.path.append("../decode_rf_data")
    sys.path.append("../tile_maps")
    from plot_tile_maps import plot_healpix
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
        default="../../outputs/tile_maps/null_test/",
        help="Output directory. Default=../../outputs/tile_maps/null_test/",
    )
    parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="../../outputs/tile_maps/tile_maps_raw/",
        help="Output directory. Default=../../outputs/tile_maps/tile_maps_raw/",
    )
    parser.add_argument(
        "--ref_model",
        metavar="\b",
        default="../../outputs/reproject_ref/ref_dipole_models.npz",
        help="Healpix reference FEE model file. default=../../outputs/reproject_ref/ref_dipole_models.npz",
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
    ref_model = args.ref_model
    nside = args.nside

    out_dir.mkdir(parents=True, exist_ok=True)

    ref_tiles = [
        "S35XX_rf0XX_sat_maps.npz",
        "S35YY_rf0YY_sat_maps.npz",
        "S35XX_rf1XX_sat_maps.npz",
        "S35YY_rf1YY_sat_maps.npz",
    ]

    za, indices = ring_indices(nside=nside)
    good_rf0XX = good_maps(ref_tiles[0])
    good_rf0YY = good_maps(ref_tiles[1])
    good_rf1XX = good_maps(ref_tiles[2])
    good_rf1YY = good_maps(ref_tiles[3])

    # Load reference FEE model
    ref_fee_model = np.load(ref_model, allow_pickle=True)
    beam_XX = ref_fee_model["XX"]
    beam_YY = ref_fee_model["YY"]

    good_rf0XX = good_rf0XX - beam_XX
    good_rf0YY = good_rf0YY - beam_YY
    good_rf1XX = good_rf1XX - beam_XX
    good_rf1YY = good_rf1YY - beam_YY

    sum_array = np.array([good_rf0XX, good_rf0YY, good_rf1XX, good_rf1YY])

    ref_median = np.nanmedian(sum_array, axis=0)

    ref_rings = []
    for i in indices:
        ref_rings.append(np.nanmedian(ref_median[i]))

    ref_med = np.median(ref_rings)

    ref_rings = np.array(ref_rings) - ref_med

    fit_ref = poly_fit(za, ref_rings, 8)

    # plt.style.use('seaborn')
    nice_fonts = {
        # Use LaTeX to write all text
        # "text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 14,
        "font.size": 14,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "ytick.color": "w",
        "xtick.color": "w",
        "axes.labelcolor": "w",
        "axes.edgecolor": "w",
    }

    plt.rcParams.update(nice_fonts)

    fig = plt.figure(figsize=(11, 6))

    plt.scatter(
        za,
        ref_rings,
        marker=".",
        s=99,
        alpha=0.9,
        color="#a7ff83",
        label=r"$\Delta$ref [$\theta$]",
    )
    plt.plot(
        za,
        fit_ref,
        linewidth=2.1,
        alpha=1,
        color="#fa4659",
        label=r"8$^{th}$ order fit",
    )

    leg = plt.legend(
        loc="lower left", frameon=True, framealpha=0.3, markerscale=1, handlelength=1
    )
    leg.get_frame().set_facecolor("white")
    for l in leg.legendHandles:
        l.set_alpha(1)
    for text in leg.get_texts():
        plt.setp(text, color="w")

    plt.xlabel(r"Zenith Angle [degrees]")
    plt.ylabel("Residual Power [dB]")
    plt.tight_layout()
    plt.savefig("ref_residuals.png", transparent=True, bbox_inches="tight")
