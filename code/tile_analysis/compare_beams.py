import sys
import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly
from mwa_pb.mwa_sweet_spots import all_grid_points

sys.path.append('../tile_maps')
from plot_tile_maps import plot_healpix
    
sys.path.append('../decode_rf_data')
from colormap import spectral, jade, kelp

# Custom spectral colormap
jade, _ = jade()

def hp_slices_horizon(nside=None):
    '''Healpix pix indices of NS, EW slices and above horizon'''
    # theta phi values of each pixel 
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    # healpix indices above the horizon
    above_horizon_indices = np.where(θ <= np.radians(85))[0]
    
    # pixel coords above the horizon
    ɸ_above_horizon = ɸ[above_horizon_indices]

    NS_indices = []
    EW_indices = []

    # pixel indices along N, E, S, W slices
    # order the indices such that they proceed from N -> S or E -> W
    n_slice = sorted(np.where((np.round(np.degrees(ɸ_above_horizon))) ==  45)[0], reverse=True)
    e_slice = sorted(np.where((np.round(np.degrees(ɸ_above_horizon))) == 135)[0], reverse=True)
    s_slice = sorted(np.where((np.round(np.degrees(ɸ_above_horizon))) == 225)[0])
    w_slice = sorted(np.where((np.round(np.degrees(ɸ_above_horizon))) == 315)[0])

    NS_indices.extend(n_slice)
    NS_indices.extend(s_slice)
    EW_indices.extend(e_slice)
    EW_indices.extend(w_slice)

    return [NS_indices, EW_indices, above_horizon_indices]


def slice_map(hp_map):
    '''slices healpix map along NS, EW'''
    
    NS_indices, EW_indices, _ = hp_slices_horizon(nside)

    θ_NS, ɸ_NS = np.degrees(hp.pix2ang(nside, NS_indices))
    θ_EW, ɸ_EW = np.degrees(hp.pix2ang(nside, EW_indices))

    zenith_angle_NS = []
    for i, j in zip(θ_NS, ɸ_NS):
        if j <= 180:
            zenith_angle_NS.append(-1*i)
        else:
            zenith_angle_NS.append(i)
    
    zenith_angle_EW = []
    for i, j in zip(θ_EW, ɸ_EW):
        if j <= 180:
            zenith_angle_EW.append(-1*i)
        else:
            zenith_angle_EW.append(i)
    
    #print(sorted(zenith_angle_NS))

    NS_data = [hp_map[NS_indices], zenith_angle_NS]
    EW_data = [hp_map[EW_indices], zenith_angle_EW]

    return [NS_data, EW_data]

def nan_mad(ref_map):
    '''Compute mad while ignoring nans'''
    ref_map_mad = []
    for j in ref_map:
        if j != []:
            j = np.asarray(j)
            j = j[~np.isnan(j)]
            ref_map_mad.append(mad(j))
        else:
            ref_map_mad.append(np.nan)

    ref_map_mad = np.asarray(ref_map_mad)
    ref_map_mad[np.where(ref_map_mad == np.nan)] = np.nanmean(ref_map_mad)

    return ref_map_mad


def map_slice(good_map):
    '''slices ref healpix map along NS & EW'''
    
    ref_map_NS, ref_map_EW = slice_map(np.asarray(good_map))

    ref_med_map_NS = np.asarray([(np.nanmedian(i) if i != [] else np.nan) for i in ref_map_NS[0]])
    # Scale mean map such that the max value is 0
    ref_med_map_scaled_NS = np.asarray([i-np.nanmax(ref_med_map_NS) for i in ref_med_map_NS])
    #ref_mad_map_NS = np.asarray([mad(i) for i in ref_map_NS[0]])
    ref_mad_map_NS = np.asarray(nan_mad(ref_map_NS[0]))
    za_NS = ref_map_NS[1]

    ref_med_map_EW = np.asarray([(np.nanmedian(i) if i != [] else np.nan) for i in ref_map_EW[0]])
    # Scale mean map such that the max value is 0
    ref_med_map_scaled_EW = np.asarray([i-np.nanmax(ref_med_map_EW) for i in ref_med_map_EW])
    #ref_mad_map_EW = np.asarray([mad(i) for i in ref_map_EW[0]])
    ref_mad_map_EW = np.asarray(nan_mad(ref_map_EW[0]))
    za_EW = ref_map_EW[1]
    
    NS_data = [ref_med_map_scaled_NS, ref_mad_map_NS, za_NS]
    EW_data = [ref_med_map_scaled_EW, ref_mad_map_EW, za_EW]

    return [NS_data, EW_data]


# rotate func written by Jack Line
def rotate(nside, angle=None,healpix_array=None,savetag=None,flip=False):
    '''Takes in a healpix array, rotates it by the desired angle, and saves it.
    Optionally flip the data, changes east-west into west-east because
    astronomy'''

    # theta phi values of each pixel 
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    new_hp_inds = hp.ang2pix(nside,θ,ɸ+angle)

    ##Flip the data to match astro conventions
    if flip == True:
        new_angles = []
        for phi in ɸ:
            if phi <= np.pi:
                new_angles.append(np.pi - phi)
            else:
                new_angles.append(3*np.pi - phi)
        new_hp_inds = hp.ang2pix(nside,ɸ, np.asarray(new_angles))

    ##Save the array in the new order
    if savetag:
        np.savez_compressed(savetag,beammap=healpix_array[new_hp_inds])

    return healpix_array[new_hp_inds]

# chisquared minimization to best fit map to data
def fit_gain(map_data=None,map_error=None,beam=None):
    '''Fit the beam model to the measured data using
    chisquared minimization'''

    bad_values = np.isnan(map_data)
    map_data = map_data[~bad_values]
    map_error = map_error[~bad_values]

    map_error[np.where(map_error == 0)] = np.mean(map_error)

    #print(len(map_data), len(map_error))
    def chisqfunc(gain):
        model = beam[~bad_values] + gain
        chisq = sum((map_data - model)**2)
        #chisq = sum(((map_data - model)/map_error)**2)
        return chisq

    x0 = np.array([0])

    result =  opt.minimize(chisqfunc,x0)
    
    return result.x


def poly_fit(x, y, map_data, order):
    '''Fit polynominal of order to data'''
    
    x = np.asarray(x)
    y = np.asarray(y)

    bad_values = np.isnan(map_data)
    x_good = x[~bad_values]
    y_good = y[~bad_values]
    coefs = poly.polyfit(x_good, y_good, order)
    fit = poly.polyval(x, coefs)
    return fit


def plt_slice(
        fig=None, sub=None,
        zen_angle=None, map_slice=None,
        map_error=None, model_slice=None, 
        delta_pow=None, pow_fit=None, 
        slice_label=None, model_label=None):

    '''Plot a slice of the beam, with measured
    data, errorbars, and fit the simulated beam
    to the data. Also plot the diff b/w data and
    the model'''

    ax = fig.add_subplot(sub)

    ax.errorbar(
           zen_angle, map_slice, yerr=map_error, 
           fmt='.', color='#326765', ecolor='#7da87b',
           elinewidth=1.2, capsize=1.2, capthick=1.4,
           alpha=0.9, label=slice_label)

    ax.plot(zen_angle,model_slice, color='#c70039', linewidth=1.2, alpha=0.9, label=model_label)

    ax.set_ylabel('Power (dB)')
    #ax.set_ylim(bottom=-30)
    #ax.set_xlabel('Zenith Angle (degrees)')
    ax.legend()
    ax.set_xlim([-90,90])
    ax.set_ylim([-32,12])

    divider = make_axes_locatable(ax)
    dax = divider.append_axes("bottom", size="30%", pad=0.1)

    #dax = fig.add_subplot(2,1,2)
    dax.scatter(zen_angle, delta_pow, marker='.', color='#27296d')
    dax.plot(zen_angle, pow_fit, linewidth=1.2, alpha=0.9, color='#ff8264')
    dax.set_ylabel('Data - Model (dB)')
    dax.set_xticklabels([])

    return ax


if __name__=='__main__':
    
    import argparse
    import numpy as np
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gs
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.stats import median_absolute_deviation as mad

    import sys
    sys.path.append('../decode_rf_data')
    from plot_tile_maps import plot_healpix
    from colormap import spectral
    
    # Custom spectral colormap
    cmap = spectral()
    
    parser = argparse.ArgumentParser(description="""
        Plot healpix map of reference data
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='../../outputs/tile_analysis/compare_beams/',
            help='Output directory. Default=../../outputs/tile_analysis/compare_beams/')
    parser.add_argument('--map_dir', metavar='\b', default='../../outputs/tile_maps/tile_maps_norm/',
            help='Tile map directory. Default=../../outputs/tile_maps/tile_maps_norm/')
    parser.add_argument('--fee_map', metavar='\b', default='../../outputs/tile_analysis/FEE_maps/mwa_fee_beam.npz',
            help='Healpix FEE map of mwa tile. default=../../outputs/tile_analysis/FEE_maps/mwa_fee_beam.npz')
    parser.add_argument('--nside', metavar='\b',type=int, default=32,help='Healpix Nside. Default = 32')
    
    
    args = parser.parse_args()
    
    out_dir     = Path(args.out_dir)
    map_dir     = Path(args.map_dir)
    fee_map     = args.fee_map
    nside       = args.nside
    
    # make output dir if it doesn't exist
    out_dir.mkdir(parents=True, exist_ok=True)

    fee_map = np.load(fee_map, allow_pickle=True)

    tile_map = np.load(f'{map_dir}/S35XX_rf0XX_tile_maps.npz', allow_pickle=True)

    pointings = ['0', '2', '4', '41']

    for p in pointings:
        fee     = fee_map[p]
        tile    = tile_map[p]
        
        # find the pointing center in radians
        pointing_center_az = np.radians(all_grid_points[int(p)][1])
        pointing_center_za = np.radians(all_grid_points[int(p)][3])
        
        # convert it to a healpix vector
        pointing_vec = hp.ang2vec(pointing_center_za, pointing_center_az)

        # find all healpix indices within 10 degrees of pointing center
        ipix_disc = hp.query_disc(nside=nside, vec=pointing_vec, radius=np.radians(10))
        
        # healpix meadian map
        tile_med = np.asarray([(np.median(i) if i != [] else np.nan ) for i in tile])
       
        # find the max value within 10 degrees of pointing center
        ipix_max = np.nanmax(tile_med[ipix_disc])
        
        # scale map such that the above max is set to 0dB 
        tile_scaled = np.asarray([(i - ipix_max) for i in tile_med])
        
        fig = plt.figure(figsize=(8,10))
        #fig.suptitle(f'Good Map: {tile}/{ref} @ {p}', fontsize=16)
        plot_healpix(data_map=tile_scaled, sub=(1,1,1), cmap=jade, vmin=-40, vmax=0)
        plt.savefig(f'test_{p}.png',bbox_inches='tight')
        plt.close()

