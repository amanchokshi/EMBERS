"""
Referece Models
---------------

Convert FEKO models on the reference antennas into usable healpix maps

"""

# Code developed by Jack Line and adapted here
# https://github.com/JLBLine/MWA_ORBCOMM

from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
import pkg_resources
from embers.rf_tools.colormaps import spectral
from embers.tile_maps.beam_utils import plot_healpix
from matplotlib import pyplot as plt
from scipy.interpolate import RectSphereBivariateSpline

matplotlib.use("Agg")
cmap, _ = spectral()


def create_model(nside, file_name=None):
    """Takes feko .ffe reference model, converts into healpix and smooths the response

    :param nside: Healpix nside
    :param file_name: Path to :samp:`.ffe` feko model of reference tiles

    :returns:
        - :class:`tuple` of (beam_response, theta_mesh, phi_mesh, power, theta)

    """

    # Make an empty array for healpix projection
    len_empty_healpix = hp.nside2npix(nside)
    beam_response = np.zeros(len_empty_healpix) * np.nan

    # Had problems with having coords = 0 when interpolating, so make a small number
    # and call it epsilon for some reason
    epsilon = 1.0e-12

    # Load ffe model data in, which is stored in real and imaginary components
    # of a phi/theta polaristaion representation of the beam
    # apparently this is normal to beam engineers

    data = np.loadtxt(file_name)
    theta = data[:, 0]
    re_theta = data[:, 2]
    im_theta = data[:, 3]
    re_phi = data[:, 4]
    im_phi = data[:, 5]

    # Convert to complex numbers
    X_theta = re_theta.astype(complex)
    X_theta.imag = im_theta

    X_phi = re_phi.astype(complex)
    X_phi.imag = im_phi

    # make a coord grid for the data in the .ffe model
    theta_range = np.linspace(epsilon, 90, 91)
    phi_range = np.linspace(epsilon, 359.0, 360)
    theta_mesh, phi_mesh = np.meshgrid(theta_range, phi_range)

    # Convert reference model into power from the complex values
    power = 10 * np.log10(abs(X_theta) ** 2 + abs(X_phi) ** 2)

    # Make an inf values super small
    power[np.where(power == -np.inf)] = -80

    # Had a problem with edge effects, leave off near horizon values
    # Basically remove the last 90 values of power list, where θ == 90
    power = power[:-91]

    # Get things into correct shape and do an interpolation
    power.shape = phi_mesh.shape

    # Bivariate spline approximation over a rectangular mesh on a sphere
    # s is a paramater I had to play with to get by eye nice results
    # s: positive smoothing factor
    lut = RectSphereBivariateSpline(
        theta_range * (np.pi / 180.0), phi_range * (np.pi / 180.0), power.T, s=0.1
    )

    # Get the theta and phi of all healpixels
    i_pix = np.arange(hp.nside2npix(nside))
    theta_rad, phi_rad = hp.pix2ang(nside, i_pix)

    # Evaluate the interpolated function at healpix gridpoints
    # Use spline to map beam into healpixels
    beam_response = lut.ev(theta_rad, phi_rad)

    # rescale beam response to have peak value of 0 dB
    beam_response = beam_response - np.nanmax(beam_response)

    return beam_response, theta_mesh, phi_mesh, power, theta


def ref_healpix_save(nside, out_dir):
    """Save and plot reference healix maps

    :param nside: Healpix nside
    :param out_dir: Path to output directory

    :returns:
        - Reference beam models have to :samp:`.npz` file in :samp:`out_dir` and plot of reference beam patterns

    """

    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Reference FEE models in XX & YY pols
    refXX = pkg_resources.resource_filename(
        "embers.kindle", "data/ref_models/MWA_reference_tile_FarField_XX.ffe"
    )
    refYY = pkg_resources.resource_filename(
        "embers.kindle", "data/ref_models/MWA_reference_tile_FarField_YY.ffe"
    )

    healpix_XX, theta_mesh, phi_mesh, power_XX, theta = create_model(
        nside, file_name=refXX
    )
    healpix_YY, theta_mesh, phi_mesh, power_YY, theta = create_model(
        nside, file_name=refYY
    )

    # Plot the things to sanity check and save results
    plt.style.use("seaborn")
    fig = plt.figure(figsize=(10, 14))

    ax3 = fig.add_subplot(2, 2, 1)
    plot_healpix(
        data_map=healpix_XX,
        sub=(2, 2, 1),
        title="Reference XX",
        vmin=-20,
        vmax=0,
        cmap=cmap,
    )
    ax3.axis("off")
    ax3.set_title("Healpix XX", x=0.4, y=-0.33)

    ax4 = fig.add_subplot(2, 2, 2)
    plot_healpix(
        data_map=healpix_YY,
        sub=(2, 2, 2),
        title="Reference YY",
        vmin=-20,
        vmax=0,
        cmap=cmap,
    )
    ax4.axis("off")
    ax4.set_title("Healpix YY", x=0.6, y=-0.33)

    np.savez_compressed(
        f"{out_dir}/ref_dipole_models.npz", XX=healpix_XX, YY=healpix_YY
    )
    fig.savefig(f"{out_dir}/reproject_dipole_models.png", bbox_inches="tight")
