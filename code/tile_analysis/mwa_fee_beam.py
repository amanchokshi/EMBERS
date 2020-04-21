import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.stats import median_absolute_deviation
from itertools import cycle

from matplotlib import rcParams
from os import environ
import scipy.optimize as opt
from mpl_toolkits.axes_grid1 import make_axes_locatable

MWAPY_H5PATH = "./../../../MWA_Beam/mwa_pb/data/mwa_full_embedded_element_pattern.h5"
from mwa_pb import primary_beam
from mwa_pb import beam_full_EE
from mwa_pb import mwa_tile
from mwa_pb.mwa_sweet_spots import all_grid_points


import sys
sys.path.append('../tile_maps')
from plot_tile_maps import plot_healpix


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
        return (np.sqrt(vis[:,:,0,0].real),np.sqrt(vis[:,:,1,1].real))
    else:
        return (vis[:,:,0,0].real,vis[:,:,1,1].real)

   


if __name__=='__main__':

    import sys
    import argparse
    from pathlib import Path
    sys.path.append('../decode_rf_data')
    from colormap import jade
    
    # Custom spectral colormap
    jade, _ = jade()

    parser = argparse.ArgumentParser(description="""
        Create Simulated MWA Beam response maps using mwapy
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/tile_analysis/FEE_maps/', help='Output directory. Default=./../../outputs/tile_analysis/FEE_maps/')
    parser.add_argument('--nside', metavar='\b', type=int,  default=32,help='Healpix Nside. Default = 32')
    
    args = parser.parse_args()
    
    out_dir     = Path(args.out_dir)
    nside       = args.nside
    
    # make output directory if it doesn't exist
    Path(out_dir).mkdir(parents=True, exist_ok=True)


    pointings = ['0', '2', '4']

    
    fee_beam = {}
    for p in pointings:

        # Empty array for model beam
        npix = hp.nside2npix(nside)
        beam_response = np.zeros(npix)

        # healpix indices above horizon
        # convert to zenith angle and azimuth
        above_horizon = range(int(npix/2))
        beam_zas,beam_azs = hp.pix2ang(nside, above_horizon)
   
        delay_p = np.array([all_grid_points[int(p)][-1], all_grid_points[int(p)][-1]])

        # S21 had a missing dipole, so need a different amplitude array for the model
        amps = np.ones((2,16))
        #if AUT == 'S21XX':
        #    amps[0,5] = 0

        # Make beam response
        response = local_beam([list(beam_zas)], [list(beam_azs)], freq=137.85e+6, delays=delay_p, zenithnorm=True, power=True, interp=False, amps=amps)
        response = response[0][0]
        
        # Stick in an array, convert to decibels, and noralise
        beam_response[above_horizon] = response
        decibel_beam = 10*np.log10(beam_response)
        normed_beam = decibel_beam - decibel_beam.max()
        fee_beam[p] = normed_beam
        
        plot_healpix(data_map=normed_beam, sub=(1,1,1), cmap=jade, vmin=-40, vmax=0)
        plt.savefig(f'{out_dir}/mwa_fee_beam_{p}.png')
        plt.close()

    np.savez_compressed(f'{out_dir}/mwa_fee_beam.npz', **fee_beam)


