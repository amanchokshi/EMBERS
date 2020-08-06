"""
Beam Utils
----------

A set of tools used to create and visualize tile maps

"""

import healpy as hp
import numpy as np


# rotate func written by Jack Line
def rotate_map(nside, angle=None, healpix_array=None, savetag=None, flip=False):
    """Rotates healpix array by the desired angle, and saves it.

    Optionally flip the data, changes east-west into west-east because astronomy

    :param nside: Healpix nside
    :param angle: Angle by which to rotate the healpix map
    :param healpix_array: Input healpix array to be rotated
    :param savetag: If given, will save the rotated beam map to :samp:`.npz` file
    :param flip: Do an astronomy coordinate flip, if True

    :returns:
        - Rotated healpix map

    """

    # theta phi values of each pixel
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    new_hp_inds = hp.ang2pix(nside, θ, ɸ + angle)

    # Flip the data to match astro conventions
    if flip is True:
        new_angles = []
        for phi in ɸ:
            if phi <= np.pi:
                new_angles.append(np.pi - phi)
            else:
                new_angles.append(3 * np.pi - phi)
        new_hp_inds = hp.ang2pix(nside, ɸ, np.asarray(new_angles))

    # Save the array in the new order
    if savetag:
        np.savez_compressed(savetag, beammap=healpix_array[new_hp_inds])

    return healpix_array[new_hp_inds]


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
    """Yeesh do some healpix magic to plot the thing

    :param data_map: Healpix input map to plot
    :param fig: Figure number to use
    :param sub: Matplotlib subplot syntax
    :param title: Plot title
    :param vmin: Colormap minimum
    :param vmax: Colormap maximum
    :param cmap: Matplotlib :class:`~matplotlib.colors.ListedColormap`
    :param cbar: If True, plot a colorbar

    :returns:
        - Plot of healpix map

    """

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
