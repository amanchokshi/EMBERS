if __name__=='__main__':
    
    import argparse
    import numpy as np
    import healpy as hp
    import matplotlib.pyplot as plt
    from pathlib import Path
    from scipy.stats import median_absolute_deviation as mad
    
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
    
    for f in out_dir.glob('*.npz'):
        f_name, _ = f.name.split('.')
        ref, _, _ = f_name.split('_')
        
        # load data from map .npz file
        map_data = np.load(f, allow_pickle=True)
        ref_map = map_data['ref_map']
        ref_counter = map_data['ref_counter']
        
        # compute the median for every pixel array
        ref_map_med = [(np.median(i) if i != [] else np.nan ) for i in ref_map]
        ref_map_mad = [mad(i) for i in ref_map]

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

        print(NS_indices)
        print(EW_indices)
        
        break

