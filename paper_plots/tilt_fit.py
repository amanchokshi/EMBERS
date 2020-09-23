from pathlib import Path

import healpy as hp
import numpy as np
from embers.tile_maps.beam_utils import plot_healpix
from healpy.rotator import euler_matrix_new as euler
from matplotlib import pyplot as plt


def reproject_map(nside, pixel, healpix_array=None):
    """Reproject beam map to be centered around healpix pixel instead of zenith.
    :param nside: Healpix nside
    :param pixel: Healpix index of new phase centre
    :param healpix_array: Input healpix array to be rotated
    :returns:
        - Reprojected healpix map
    """

    # coordinated of new phase centre
    #  theta, phi = hp.pix2ang(nside, pixel)
    # Theta moves primary beam eastward
    # Phi rotates primary beam anticlockwise. Or, E --> N
    #  theta, phi = [1, -45]
    theta, phi = hp.pix2ang(nside, pixel)

    vec = hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))
    eu_mat = euler(-phi, -theta, phi, deg=False)
    rot_map = hp.rotator.rotateVector(eu_mat, vec)
    new_hp_inds = hp.vec2pix(nside, rot_map[0], rot_map[1], rot_map[2])

    return healpix_array[new_hp_inds]


def pointing_residual(nside, pixel, beam=None, fee=None):

    # Select primary and secondary lobes
    mask = np.where(fee < -30)

    rot_beam = reproject_map(nside, pixel, healpix_array=beam)

    resi = rot_beam - fee
    resi[mask] = np.nan

    # Select central 40 deg of beam. Primary + edges of secondary lobes
    ind = int(hp.ang2pix(nside, np.radians(40), 0))
    resi = resi[:ind]

    # Mean of residuals as proxy for pointing errors
    resi_mean = np.nanmean(resi)

    return(resi_mean)

if __name__ == "__main__":
    pass
