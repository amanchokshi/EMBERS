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


if __name__=='__main__':
    
    import argparse
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    from scipy.stats import median_absolute_deviation as mad

    import matplotlib.pyplot as plt
    
    import sys
    sys.path.append('../decode_rf_data')
    from colormap import spectral
    
    # Custom spectral colormap
    cmap = spectral()
    
    parser = argparse.ArgumentParser(description="""
        Plot healpix map of reference data
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='../../outputs/null_test/',help='Output directory. Default=../../outputs/null_test/')
    parser.add_argument('--nside', metavar='\b', default=32,help='Healpix Nside. Default = 32')
    
    
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    nside   = args.nside
    
    ref_tiles = ['rf0XX', 'rf0YY', 'rf1XX', 'rf1YY']

    rf0XX_NS, rf0XX_EW = ref_map_slice(ref_tiles[0])
    rf0YY_NS, rf0YY_EW = ref_map_slice(ref_tiles[1])
    rf1XX_NS, rf1XX_EW = ref_map_slice(ref_tiles[2])
    rf1YY_NS, rf1YY_EW = ref_map_slice(ref_tiles[3])
    

