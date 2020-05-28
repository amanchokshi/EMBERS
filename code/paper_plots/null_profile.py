import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly
from scipy.stats import median_absolute_deviation as mad


def ring_indices(nside=None):
    '''Healpix pix indices of NS, EW slices and above horizon'''
    # theta phi values of each pixel 
    hp_indices = np.arange(hp.nside2npix(nside))
    θ, ɸ = hp.pix2ang(nside, hp_indices)

    za = np.arange(0,87,2)

    ring_indices = []
    for i in range(len(za)):
        if i < 43:
            r_indices = np.where(np.logical_and(θ >= np.radians(za[i]), θ < np.radians(za[i+1])))[0]
            ring_indices.append(r_indices)
    za_cen = za[:-1] + 1
    return [za_cen, ring_indices]

def good_maps(ref_map):
    '''Creates a ref map with only good satellites'''
   
    pointings = ['0','2','4']
    
    # load data from map .npz file
    f = Path(f'{map_dir}/{ref_map}')
    tile_data = np.load(f, allow_pickle=True)
    tile_data = {key:tile_data[key].item() for key in tile_data}
    ref_map = tile_data['ref_map']  
    
    # Good sats from which to make plots
    good_sats = [
            25338, 25982, 25984, 25985,
            28654, 40086, 40087,
            40091, 41179, 41180,
            41182, 41183, 41184,
            41185, 41187, 41188,
            41189, 44387
            ]
    
    orbcomm = [
            25982, 25984, 25985,
            40086, 40087, 40091,
            41179, 41180, 41182, 
            41183, 41184, 41185, 
            41187, 41188, 41189
            ]
    noaa = [25338, 28654]

    meteor = [44387]
    
    # Empty good map
    good_map = [[] for pixel in range(hp.nside2npix(nside))]
    
    for p in pointings:
        
        # append to good map from all good sat data
        for sat in meteor:
            for pix in range(hp.nside2npix(nside)):
                good_map[pix].extend(ref_map[p][sat][pix])
        
    mad_map = []
    for j in good_map:
        if j != []:
            j = np.asarray(j)
            j = j[~np.isnan(j)]
            mad_map.append(mad(j))
        else:
            mad_map.append(np.nan)
    
    good_map = [np.nanmedian(pixel) for pixel in good_map]
        
    return (good_map, mad_map)

def poly_fit(x, y, order):
    '''Fit polynominal of order to data'''
    
    x = np.asarray(x)
    y = np.asarray(y)

    coefs = poly.polyfit(x, y, order)
    fit = poly.polyval(x, coefs)
    return fit


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
    sys.path.append('../tile_maps')
    from plot_tile_maps import plot_healpix
    from colormap import spectral
    
    # Custom spectral colormap
    cmap = spectral()
    
    parser = argparse.ArgumentParser(description="""
        Plot healpix map of reference data
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='../../outputs/tile_maps/null_test/',
            help='Output directory. Default=../../outputs/tile_maps/null_test/')
    parser.add_argument('--map_dir', metavar='\b', default='../../outputs/tile_maps/tile_maps_raw/',
            help='Output directory. Default=../../outputs/tile_maps/tile_maps_raw/')
    parser.add_argument('--ref_model', metavar='\b', default='../../outputs/reproject_ref/ref_dipole_models.npz',
            help='Healpix reference FEE model file. default=../../outputs/reproject_ref/ref_dipole_models.npz')
    parser.add_argument('--nside', metavar='\b',type=int, default=32,help='Healpix Nside. Default = 32')
    
    
    args = parser.parse_args()
    
    out_dir     = Path(args.out_dir)
    map_dir     = Path(args.map_dir)
    ref_model   = args.ref_model
    nside       = args.nside
    
    out_dir.mkdir(parents=True, exist_ok=True)
    
    ref_tiles = [
            'S35XX_rf0XX_sat_maps.npz', 
            'S35YY_rf0YY_sat_maps.npz', 
            'S35XX_rf1XX_sat_maps.npz', 
            'S35YY_rf1YY_sat_maps.npz'
            ]

    
    za, indices = ring_indices(nside=nside)
    
    # median and mad values of induvidual ref maps 
    good_rf0XX, mad_rf0XX = good_maps(ref_tiles[0])
    good_rf0YY, mad_rf0YY = good_maps(ref_tiles[1])
    good_rf1XX, mad_rf1XX = good_maps(ref_tiles[2])
    good_rf1YY, mad_rf1YY = good_maps(ref_tiles[3])

    # Load reference FEE model
    ref_fee_model = np.load(ref_model, allow_pickle=True)
    beam_XX = ref_fee_model['XX']
    beam_YY = ref_fee_model['YY']
   
    # residuals of individual maps
    res_rf0XX = good_rf0XX - beam_XX
    res_rf0YY = good_rf0YY - beam_YY
    res_rf1XX = good_rf1XX - beam_XX
    res_rf1YY = good_rf1YY - beam_YY

    sum_array = np.array([res_rf0XX, res_rf0YY, res_rf1XX, res_rf1YY])
    mad_array = np.array([mad_rf0XX, mad_rf0YY, mad_rf1XX, mad_rf1YY])
    
    # sum along corresponding pixels
    res_median = np.nanmedian(sum_array, axis=0)
    mad_average = np.nanmean(mad_array, axis=0)
   
    # median of values in zenith angle rings
    res_rings = []
    mad_rings = []
    for i in indices:
        res_rings.append(np.nanmedian(res_median[i]))
        mad_rings.append(np.nanmedian(mad_average[i]))

    # Scale the residuals to 0 median 
    res_med = np.median(res_rings)
    res_rings = np.array(res_rings) - res_med

    fit_res = poly_fit(za, res_rings, 8)  
    
    plt.style.use('seaborn')
    nice_fonts = {
            # Use LaTeX to write all text
            #"text.usetex": True,
            "font.family": "sans-serif",
            # Use 10pt font in plots, to match 10pt font in document
            "axes.labelsize": 8,
            "font.size": 10,
            # Make the legend/label fonts a little smaller
            "legend.fontsize": 6,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            }
    
    plt.rcParams.update(nice_fonts)
    
    fig = plt.figure(figsize=(3.6,2.4))
           
    plt.style.use('seaborn')
    plt.errorbar(
           za, res_rings, yerr=mad_rings, 
           fmt='.', color='#326765', ecolor='#7da87b',
           elinewidth=1.4, capsize=1.4, capthick=1.6,
           alpha=0.9, ms=9, label=r'$\Delta$ref [$\theta$]')
    
    #plt.scatter(za, res_rings, marker='.', s=77, alpha=0.9,  color='#466844', label=r'$\Delta$ref [$\theta$]')
    plt.plot(za, fit_res, linewidth=1.4, alpha=1, color='#fa4659', label=r'8$^{th}$ order fit')
    
    leg = plt.legend(frameon=True, markerscale=1, handlelength=1)
    leg.get_frame().set_facecolor('white')
    for l in leg.legendHandles:
        l.set_alpha(1)
    
    plt.xlabel(r'Zenith Angle [degrees]')
    plt.ylabel('Residual Power [dB]')
    plt.tight_layout()
    #plt.savefig('../../outputs/paper_plots/ref_residuals.pdf', bbox_inches='tight')
    plt.savefig('ref_residuals_meteor.pdf', bbox_inches='tight')



