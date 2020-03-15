import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly


def hp_slices_horizon(nside=None):
    '''Healpix pix indices of NS, EW slices and above horizon'''
    # theta phi values of each pixel 
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    # healpix indices above the horizon
    above_horizon_indices = np.where(θ <= np.radians(90))[0]
    
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
    θ_EW, ɸ_EW = hp.pix2ang(nside, EW_indices)

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


def ref_map_slice(ref_tile):
    '''slices ref healpix map along NS, EW'''
    
    # load data from map .npz file
    map_data = np.load(f'{out_dir}/{ref_tile}_map_healpix.npz', allow_pickle=True)
    ref_map = map_data['ref_map']
    ref_counter = map_data['ref_counter']
    
    # compute the median for every pixel array
    ref_map_med = np.asarray([(np.median(i) if i != [] else np.nan ) for i in ref_map])
    ref_map_mad = np.asarray([mad(i) for i in ref_map])

    # slice the median map along NS, EW
    ref_med_NS, ref_med_EW = slice_map(ref_map_med)
    ref_med_map_NS, za_NS = ref_med_NS
    ref_med_map_EW, za_EW = ref_med_EW

    # slice the mad map along NS, EW
    ref_mad_NS, ref_mad_EW = slice_map(ref_map_mad)
    ref_mad_map_NS, _ = ref_mad_NS
    ref_mad_map_EW, _ = ref_mad_EW
    
    NS_data = [ref_med_map_NS, ref_mad_map_NS, za_NS]
    EW_data = [ref_med_map_EW, ref_mad_map_EW, za_EW]

    return [NS_data, EW_data]


# rotate func written by Jack Line
def rotate(angle=None,healpix_array=None,savetag=None,flip=False):
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
        savez_compressed(savetag,beammap=healpix_array[new_hp_inds])

    return healpix_array[new_hp_inds]

# chisquared minimization to best fit map to data
def fit_gain(map_data=None,map_error=None,beam=None):

    bad_values = np.isnan(map_data)
    map_data = map_data[~bad_values]
    map_error = map_error[~bad_values]

    map_error[np.where(map_error == 0)] = np.mean(map_error)

    def chisqfunc(gain):
        model = beam[~bad_values] + gain
        chisq = sum(((map_data - model)/map_error)**2)
        return chisq

    x0 = np.array([0])

    result =  opt.minimize(chisqfunc,x0)
    #print(f'chisq is:{chisqfunc(result.x)}')
    #print(f'gain is: {result.x}')
    #print result.x
    return result.x


def poly_fit(x, y, order):
    coefs = poly.polyfit(x, y, order)
    fit = poly.polyval(x, coefs)

    return fit


def plt_slice(sub=None, title=None, zen_angle=None, map_slice=None, map_error=None, model_slice=None, delta_pow=None, pow_fit=None, slice_label=None, model_label=None):

    ax = fig.add_subplot(sub)

    ax.errorbar(
           zen_angle, map_slice, yerr=map_error, 
           fmt='.', color='#326765', ecolor='#7da87b',
           elinewidth=1.2, capsize=1.2, capthick=1.4,
           alpha=0.9, label=slice_label)

    ax.plot(zen_angle,model_slice, color='#c70039', linewidth=1.2, alpha=0.9, label=model_label)

    ax.set_ylabel('Power (dB)')
    ax.set_xlabel('Zenith Angle (degrees)')
    ax.legend()

    divider = make_axes_locatable(ax)
    dax = divider.append_axes("bottom", size="30%", pad=0.1)

    #dax = fig.add_subplot(2,1,2)
    dax.scatter(zen_angle, delta_pow, marker='.', color='#27296d')
    dax.plot(zen_angle, pow_fit, linewidth=1.2, alpha=0.9, color='#ff8264')
    dax.set_ylabel('Data - Model (dB)')

    return ax



if __name__=='__main__':
    
    import argparse
    import numpy as np
    from pathlib import Path
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.stats import median_absolute_deviation as mad

    import sys
    sys.path.append('../decode_rf_data')
    from plot_healpix import plot_healpix
    from colormap import spectral
    
    # Custom spectral colormap
    cmap = spectral()
    
    parser = argparse.ArgumentParser(description="""
        Plot healpix map of reference data
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='../../outputs/null_test/',
            help='Output directory. Default=../../outputs/null_test/')
    parser.add_argument('--ref_model', metavar='\b', default='../../outputs/reproject_ref/ref_dipole_models.npz',
            help='Healpix reference FEE model file. default=../../outputs/reproject_ref/ref_dipole_models.npz')
    parser.add_argument('--nside', metavar='\b', default=32,help='Healpix Nside. Default = 32')
    
    
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    ref_model = args.ref_model
    nside   = args.nside
    
    ref_tiles = ['rf0XX', 'rf0YY', 'rf1XX', 'rf1YY']

    rf0XX_NS, rf0XX_EW = ref_map_slice(ref_tiles[0])
    rf0YY_NS, rf0YY_EW = ref_map_slice(ref_tiles[1])
    rf1XX_NS, rf1XX_EW = ref_map_slice(ref_tiles[2])
    rf1YY_NS, rf1YY_EW = ref_map_slice(ref_tiles[3])

    ref_fee_model = np.load(ref_model, allow_pickle=True)
    beam_XX = ref_fee_model['XX']
    beam_YY = ref_fee_model['YY']
   
    # Rotate beam models by pi/4 to match the rotated ref data
    rotated_XX = rotate(angle=+(1*np.pi)/4.0,healpix_array=beam_XX)
    rotated_YY = rotate(angle=+(1*np.pi)/4.0,healpix_array=beam_YY)
    
    # These plots show that the pi/4 rotation was correct
    #plot_healpix(data_map=rotated_XX,sub=(1,1,1), cmap=cmap, vmin=-40, vmax=-20)
    #plot_healpix(data_map=rotated_YY,sub=(1,1,1), cmap=cmap, vmin=-40, vmax=-20)
    #plt.show()
    
    # slice the XX rotated map along NS, EW
    XX_NS, XX_EW = slice_map(rotated_XX)
    XX_NS_slice, za_NS = XX_NS
    XX_EW_slice, za_EW = XX_EW

    # slice the YY rotated map along NS, EW
    YY_NS, YY_EW = slice_map(rotated_YY)
    YY_NS_slice, za_NS = YY_NS
    YY_EW_slice, za_EW = YY_EW

    gain = fit_gain(rf0XX_NS[0], rf0XX_NS[1], XX_NS_slice)
    XX_NS_slice_norm = XX_NS_slice + gain

    del_p = rf0XX_NS[0] - XX_NS_slice_norm

    fit = poly_fit(za_NS, del_p, 3)


    plt.style.use('seaborn')
    fig = plt.figure(figsize=(8,10))

    ax1 = plt_slice(
            sub=221, title='rf0XX NS Slice', 
            zen_angle=za_NS, map_slice=rf0XX_NS[0],
            map_error=rf0XX_NS[1], model_slice=XX_NS_slice_norm,
            delta_pow=del_p, pow_fit=fit,
            slice_label='rf0XX NS', model_label='FEE Model')


    plt.tight_layout()
    plt.show()

