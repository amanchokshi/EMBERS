import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.stats import median_absolute_deviation
from itertools import cycle

from matplotlib import rcParams
from os import environ
import scipy.optimize as opt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from sys import path
path.append('./../../../MWA_Beam/mwa_pb/mwa_pb')
import primary_beam
import beam_full_EE
import mwa_tile
MWAPY_H5PATH = "./../../../MWA_Beam/mwa_pb/mwa_pb/data/mwa_full_embedded_element_pattern.h5"


def local_beam(za, az, freq, delays=None, zenithnorm=True, power=True, jones=False, interp=True, pixels_per_deg=5,amps=None):
    '''Code pulled from mwapy that generates the MWA beam response - removes unecessary extra code from mwapy/pb
        - delays is a 2x16 array, with the first 16 delays for the XX, 
          second 16 for the YY pol. values match what you find in the metafits file
        - amps are the amplitudes of each individual dipole, again in a 2x16,
            with XX first then YY'''
    tile=beam_full_EE.ApertureArray(MWAPY_H5PATH,freq)
    mybeam=beam_full_EE.Beam(tile, delays,amps)
    if interp:
        j=mybeam.get_interp_response(az, za, pixels_per_deg)
    else:
        j=mybeam.get_response(az, za)
    if zenithnorm==True:
        j=tile.apply_zenith_norm_Jones(j) #Normalise

    #Use swapaxis to place jones matrices in last 2 dimensions
    #insead of first 2 dims.
    if len(j.shape)==4:
        j=np.swapaxes(np.swapaxes(j,0,2),1,3)
    elif len(j.shape)==3: #1-D
        j=np.swapaxes(np.swapaxes(j,1,2),0,1)
    else: #single value
        pass

    if jones:
        return j

    #Use mwa_tile makeUnpolInstrumentalResponse because we have swapped axes
    vis = mwa_tile.makeUnpolInstrumentalResponse(j,j)
    if not power:
        return (sqrt(vis[:,:,0,0].real),sqrt(vis[:,:,1,1].real))
    else:
        return (vis[:,:,0,0].real,vis[:,:,1,1].real)

   


if __name__=='__main__':

    import argparse
    from pathlib import Path

    parser = argparse.ArgumentParser(description="""
        Create Simulated MWA Beam response maps using mwapy
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/tile_analysis/FEE_maps/', help='Output directory. Default=./../../outputs/tile_analysis/FEE_maps/')
    parser.add_argument('--nside', metavar='\b', type=int,  default=32,help='Healpix Nside. Default = 32')
    
    args = parser.parse_args()
    
    out_dir     = Path(args.outdir)
    nside       = args.nside

    # Empty array for model beam
    beam_response = np.zeros(hp.nside2npix(nside))
    
    beam_zas,beam_azs = hp.pix2ang(nside, all_indexes)
    beam_azs[beam_azs < 0] += 2*np.pi
    beam_azs -= np.pi / 4.0
    
    # S21 had a missing dipole, so need a different amplitude array for the model
    amps = np.ones((2,16))
    #if AUT == 'S21XX':
    #    amps[0,5] = 0

    # Make beam response
    response = local_beam([list(beam_zas)], [list(beam_azs)], freq=137.85e+6, delays=np.zeros((2,16)), zenithnorm=True, power=True, interp=False, amps=amps)
    response = response[0][0]
    
    # Stick in an array, convert to decibels, and noralise
    beam_response[all_indexes] = response
    decibel_beam = 10*np.log10(beam_response)
    normed_beam = decibel_beam - decibel_beam.max()
    
    # Just plot it above 20 deg elevation
    above_thresh = np.where(beam_zas <= 20.0*(np.pi/180.0))

