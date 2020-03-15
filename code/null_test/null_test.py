import numpy as np
import healpy as hp


def hp_slices_horizon(nside=None):
    '''Healpix pix indices of NS, EW slices and above horizon'''
    # theta phi values of each pixel 
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    # healpix indices above the horizon
    above_horizon_indices = np.where(θ <= np.radians(90))[0]
    
    # pixel coords above the horizon
    θ_above_horizon = θ[above_horizon_indices]
    ɸ_above_horizon = ɸ[above_horizon_indices]

    NS_indices = []
    EW_indices = []

    # pixel indices along N, E, S, W slices
    n_slice = np.where((np.round(np.degrees(ɸ_above_horizon))) ==  45)[0]
    e_slice = np.where((np.round(np.degrees(ɸ_above_horizon))) == 135)[0]
    s_slice = np.where((np.round(np.degrees(ɸ_above_horizon))) == 225)[0]
    w_slice = np.where((np.round(np.degrees(ɸ_above_horizon))) == 315)[0]

    NS_indices.extend(n_slice)
    NS_indices.extend(s_slice)
    EW_indices.extend(e_slice)
    EW_indices.extend(w_slice)
    
    NS_indices = sorted(NS_indices)
    EW_indices = sorted(EW_indices)
    
    return [NS_indices, EW_indices, above_horizon_indices]


def ref_map_slice(ref_tile):
    '''slices healpix map along NS, EW'''
    
    # load data from map .npz file
    map_data = np.load(f'{out_dir}/{ref_tile}_map_healpix.npz', allow_pickle=True)
    ref_map = map_data['ref_map']
    ref_counter = map_data['ref_counter']
    
    # compute the median for every pixel array
    ref_map_med = np.asarray([(np.median(i) if i != [] else np.nan ) for i in ref_map])
    ref_map_mad = np.asarray([mad(i) for i in ref_map])

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


    NS_data = [ref_map_med[NS_indices], ref_map_mad[NS_indices], zenith_angle_NS]
    EW_data = [ref_map_med[EW_indices], ref_map_mad[EW_indices], zenith_angle_EW]

    return [NS_data, EW_data]


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





if __name__=='__main__':
    
    import argparse
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    from scipy.stats import median_absolute_deviation as mad

    import matplotlib.pyplot as plt
    
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
    
    rotate_beam_XX = rotate(angle=+(1*np.pi)/4.0,healpix_array=beam_XX)
    
    plot_healpix(data_map=rotate_beam_XX,sub=(1,1,1), cmap=cmap, vmin=-40, vmax=-20)
    plt.show()


    ## Try a number of rotations to get the reference models into the same frame
    ## as the data. Feel free to double check all this
    #rotate_ungeastern = rotate(angle=+(3*pi)/4.0,healpix_array=ung_eastern,savetag='rotated+180_ung_eastern.npz')
    #rotate_ungwestern = rotate(angle=+(3*pi)/4.0,healpix_array=ung_western,savetag='rotated+180_ung_western.npz')
    #
    #rotate_ungeastern = rotate(angle=+pi/4.0,healpix_array=ung_eastern,savetag='rotated_ung_eastern.npz')
    #rotate_ungwestern = rotate(angle=+pi/4.0,healpix_array=ung_western,savetag='rotated_ung_western.npz')
    #
    #rotate_ungeastern = rotate(angle=+pi/4.0,healpix_array=ung_eastern,savetag='rotated+flip_ung_eastern.npz',flip=True)
    #rotate_ungwestern = rotate(angle=+pi/4.0,healpix_array=ung_western,savetag='rotated+flip_ung_western.npz',flip=True)
    #
    #rotate_rf0 = rotate(angle=0,healpix_array=med_dB_rf0,savetag='prerotated_rf0.npz')
    #rotate_rf1 = rotate(angle=0,healpix_array=med_dB_rf1,savetag='prerotated_rf1.npz')
    #
    #rotate_rf0_mad = rotate(angle=0,healpix_array=mad_dB_rf0,savetag='prerotated_rf0_error.npz')
    #rotate_rf1_mad = rotate(angle=0,healpix_array=mad_dB_rf1,savetag='prerotated_rf1_error.npz')


    #plt.style.use('seaborn')
    #plt.errorbar(
    #        rf0XX_NS[2], rf0XX_NS[0], yerr=rf0XX_NS[1], 
    #        fmt='o', color='#326765', ecolor='#7da87b',
    #        elinewidth=1.2, capsize=2, capthick=1.4,
    #        alpha=0.9, label='rf0XX NS')
    #plt.ylabel('Power (dB)')
    #plt.xlabel('Zenith Angle (degrees)')
    #plt.legend()
    ##plt.tight_layout()
    #plt.show()
    

