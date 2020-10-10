import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import jade, spectral
from embers.tile_maps.beam_utils import (healpix_cardinal_slices, plot_healpix,
                                         rotate_map, poly_fit)
from healpy.rotator import euler_matrix_new as euler
from matplotlib import pyplot as plt


def reproject_map(nside, phi, healpix_array=None):
    """Reproject beam map to be centered around healpix pixel instead of zenith.

    :param nside: Healpix nside
    :param phi: Angle to rotate beam by, from E --> N (anticlockwise)
    :param healpix_array: Input healpix array to be rotated
    :returns:
        - Reprojected healpix map
    """

    vec = hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))
    eu_mat = euler(-phi, 0, 0, deg=True)
    rot_map = hp.rotator.rotateVector(eu_mat, vec)
    new_hp_inds = hp.vec2pix(nside, rot_map[0], rot_map[1], rot_map[2])

    return healpix_array[new_hp_inds]


if __name__ == "__main__":

    spec, _ = spectral()
    jade, _ = jade()
    nside = 512

    ref_model = "ref_rot/ref_dipole_models.npz"

    # Load reference FEE model
    # Rotate the fee models by -pi/2 to move model from spherical (E=0) to Alt/Az (N=0)
    ref_fee_model = np.load(ref_model, allow_pickle=True)

    ref_xx = rotated_fee = rotate_map(
        nside, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee_model["XX"]
    )
    ref_yy = rotated_fee = rotate_map(
        nside, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee_model["YY"]
    )

    ref_xx_rot = reproject_map(nside, 2, healpix_array=ref_xx)

    res = ref_xx - ref_xx_rot

    NS, EW = healpix_cardinal_slices(nside, res, 85)

    # Plotting stuff
    plt.style.use("seaborn")

    nice_fonts = {
        "font.family": "sans-serif",
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "font.size": 8,
        "legend.fontsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }

    plt.rcParams.update(nice_fonts)

    plt.figure(figsize=(6, 5))

    plt.plot(NS[1], NS[0], label="NS_XX")
    plt.plot(EW[1], EW[0], label="EW_XX")
    plt.plot(NS[1][::36], poly_fit(NS[1][::36], NS[0][::36], NS[0][::36], 9), label="NS_XX")
    plt.plot(EW[1][::36], poly_fit(EW[1][::36], EW[0][::36], EW[0][::36], 9), label="EW_XX")
    plt.xlabel("Zenith Angle [degrees]")
    plt.ylabel("Residual Power [dB]")

    leg = plt.legend(loc="lower center", frameon=True, markerscale=4, handlelength=1)
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.savefig("ref_rot/res_slices_2.pdf")
    plt.close()

    #  plot_healpix(data_map=res, cmap=spec)
    #  plt.savefig("ref_rot/ref_res_2.png")
