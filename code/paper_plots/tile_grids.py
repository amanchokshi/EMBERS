import sys
import numpy as np
import healpy as hp
import scipy.optimize as opt
import numpy.polynomial.polynomial as poly
from mwa_pb.mwa_sweet_spots import all_grid_points

sys.path.append('../decode_rf_data')
from colormap import spectral, jade, kelp

# Custom spectral colormap
jade, _ = jade()


def plot_healpix(data_map=None, fig=None, sub=None,title=None,vmin=None,vmax=None,cmap=None, cbar=True):
    '''Yeesh do some healpix magic to plot the thing'''
    
    # Disable cryptic healpy warnings. Can't figure out where they originate
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    
    hp.delgraticules()
    hp.orthview(
            map=data_map,coord='E', fig=fig,
            half_sky=True, rot=(0,90,180), xsize=1200,
            title=title,sub=sub,min=vmin,max=vmax,
            cmap=cmap,notext=True, hold=True, cbar=cbar, return_projected_map=False)

    hp.graticule(dpar=10,coord='E',color='k',alpha=0.7,dmer=45, lw=0.4, ls=':')
   
    # Altitude grid
    hp.projtext(00.0*(np.pi/180.0), 225.0*(np.pi/180),'0',  color='k', coord='E', fontsize=6, fontweight='light')
    hp.projtext(30.0*(np.pi/180.0), 225.0*(np.pi/180),'30', color='k', coord='E', fontsize=6, fontweight='light')
    hp.projtext(60.0*(np.pi/180.0), 225.0*(np.pi/180),'60', color='k', coord='E', fontsize=6, fontweight='light')

    # NSEW 
    hp.projtext(80.0*(np.pi/180.0), 000.0*(np.pi/180.0), r'$N  $', coord='E',color='w', verticalalignment='top')
    hp.projtext(80.0*(np.pi/180.0), 090.0*(np.pi/180.0), r'$E  $', coord='E',color='w', horizontalalignment='right')
    hp.projtext(80.0*(np.pi/180.0), 180.0*(np.pi/180.0), r'$S  $', coord='E',color='w', verticalalignment='bottom')
    hp.projtext(80.0*(np.pi/180.0), 270.0*(np.pi/180.0), r'$W  $', coord='E',color='w', horizontalalignment='left')


def beam_maps(f):
    '''Returns pairs of [AUT,FEE] maps, for each pointing'''
    f_name, _ = Path(f).name.split('.')
    t_name, r_name, _, _ = f_name.split('_')

    #pointings = ['0','2','4','41']
    pointings = ['0','2','4']
    
    maps = []

    # load data from map .npz file
    tile_map = np.load(f, allow_pickle=True)
    fee_m = np.load(fee_map, allow_pickle=True)
    
    for p in pointings:

        tile = tile_map[p]

        if 'XX' in t_name:
            fee = fee_m[p][0]
        else:
            fee = fee_m[p][1]
        

        # Visualize the tile map and diff map
        # healpix meadian map
        tile_med = np.asarray([(np.nanmedian(j) if j != [] else np.nan ) for j in tile])
       
        maps.append([tile_med, fee])
        
        #residuals = tile_med - fee
        #residuals[np.where(fee < -30)] = np.nan
        #residuals[np.where(tile_med == np.nan)] = np.nan

    return maps


def plt_grid(f_xx, f_yy):
    # This is an Awesome plot

    maps_xx = beam_maps(f_xx) 
    maps_yy = beam_maps(f_yy)
    
    nice_fonts = {
            # Use LaTeX to write all text
            #"text.usetex": True,
            "font.family": "sans-serif",
            # Use 10pt font in plots, to match 10pt font in document
            "axes.labelsize": 8,
            "font.size": 8,
            # Make the legend/label fonts a little smaller
            "legend.fontsize": 6,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            }
    plt.rcParams.update(nice_fonts)
    #plt.rcParams['grid.color'] = '#cccccc'
    #plt.rcParams['grid.linestyle'] = ':'
    #plt.rcParams['grid.linewidth'] = 0.3
    #plt.style.use('seaborn')
   
    fig1 = plt.figure(figsize=(7.6,9.0))

    ax1  = fig1.add_axes([0.60, 0.75, 0.29, 0.22])
    plot_healpix(data_map=maps_xx[2][0],  sub=(4,3,1),  fig=fig1, title='S07XX @ 4', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax1 = plt.gca()
    image = ax1.get_images()[0]
    cax1 = fig1.add_axes([0.91, 0.75, 0.015, 0.22])
    cbar1 = fig1.colorbar(image, cax=cax1, label='Power [dB]')
    
    ax2  = fig1.add_axes([0.30, 0.75, 0.29, 0.22])
    plot_healpix(data_map=maps_xx[1][0],  sub=(4,3,2),  fig=fig1, title='S07XX @ 2', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax3  = fig1.add_axes([0.01, 0.75, 0.29, 0.22])
    plot_healpix(data_map=maps_xx[0][0],  sub=(4,3,3),  fig=fig1, title='S07XX @ Zenith', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax4  = fig1.add_axes([0.60, 0.50, 0.29, 0.22])
    plot_healpix(data_map=maps_xx[2][1],  sub=(4,3,4),  fig=fig1, title='FEE XX @ 4', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax4 = plt.gca()
    image = ax4.get_images()[0]
    cax2 = fig1.add_axes([0.91, 0.50, 0.015, 0.22])
    cbar2 = fig1.colorbar(image, cax=cax2, label='Power [dB]')
    
    ax5  = fig1.add_axes([0.30, 0.50, 0.29, 0.22])
    plot_healpix(data_map=maps_xx[1][1],  sub=(4,3,5),  fig=fig1, title='FEE XX @ 2', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax6  = fig1.add_axes([0.01, 0.50, 0.29, 0.22])
    plot_healpix(data_map=maps_xx[0][1],  sub=(4,3,6),  fig=fig1, title='FEE XX @ Zenith', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax7  = fig1.add_axes([0.60, 0.25, 0.29, 0.22])
    plot_healpix(data_map=maps_yy[2][0], sub=(4,3,7),  fig=fig1, title='S07YY @ 4', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax7 = plt.gca()
    image = ax7.get_images()[0]
    cax3 = fig1.add_axes([0.91, 0.25, 0.015, 0.22])
    cbar3 = fig1.colorbar(image, cax=cax3, label='Power [dB]')
    
    ax8  = fig1.add_axes([0.30, 0.25, 0.29, 0.22])
    plot_healpix(data_map=maps_yy[1][0], sub=(4,3,8),  fig=fig1, title='S07YY @ 2', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax9  = fig1.add_axes([0.01, 0.25, 0.29, 0.22])
    plot_healpix(data_map=maps_yy[0][0],  sub=(4,3,9),  fig=fig1, title='S07YY @ Zenith', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax10 = fig1.add_axes([0.60, 0.0, 0.29, 0.22])
    plot_healpix(data_map=maps_yy[2][1],  sub=(4,3,10), fig=fig1, title='FEE YY @ 4', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax10 = plt.gca()
    image = ax10.get_images()[0]
    cax4 = fig1.add_axes([0.91, 0.0, 0.015, 0.22])
    cbar4 = fig1.colorbar(image, cax=cax4, label='Power [dB]')
    
    ax11 = fig1.add_axes([0.30, 0.0, 0.29, 0.22])
    plot_healpix(data_map=maps_yy[1][1], sub=(4,3,11), fig=fig1, title='FEE YY @ 2', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    ax12 = fig1.add_axes([0.01, 0.0, 0.29, 0.22])
    plot_healpix(data_map=maps_yy[0][1], sub=(4,3,12), fig=fig1, title='FEE YY @ Zenith', cmap=jade, vmin=-50, vmax=0,  cbar=False)
    
    plt.savefig('test.pdf', bbox_inches='tight', dpi=300)


#    fig1 = plt.figure(figsize=(10, 8))
#
#    ax1 = plt.subplot(2,2,1)
#    plt_slice(
#            fig=fig1, sub=221,
#            zen_angle=NS_tile[2], map_slice=NS_tile_med,
#            map_error=NS_tile[1], model_slice=NS_fee[0],
#            delta_pow=del_NS, pow_fit=fit_NS, 
#            slice_label='Tile NS', model_label='FEE NS')
#
#    ax2 = fig1.add_axes([0.48, 0.52, 0.48, 0.43])
#    plot_healpix(data_map=tile_med, sub=(2,2,2), fig=fig1, title='tile map', cmap=jade, vmin=-50, vmax=0, cbar=False)
#    ax7 = plt.gca()
#    image = ax7.get_images()[0]
#    cax = fig1.add_axes([0.92, 0.52, 0.015, 0.43])
#    cbar = fig1.colorbar(image, cax=cax, label='dB')
#    
#    ax3 = plt.subplot(2,2,3)
#    plt_slice(
#            fig=fig1, sub=223,
#            zen_angle=EW_tile[2], map_slice=EW_tile_med,
#            map_error=EW_tile[1], model_slice=EW_fee[0],
#            delta_pow=del_EW, pow_fit=fit_EW, 
#            slice_label='Tile EW', model_label='FEE EW')
#    
#    ax4 = fig1.add_axes([0.48, 0.02, 0.48, 0.43])
#    plot_healpix(data_map=residuals, sub=(2,2,4), fig=fig1, title='diff map', cmap='inferno', vmin=-10, vmax=5,  cbar=False)
#    ax8 = plt.gca()
#    image = ax8.get_images()[0]
#    cax = fig1.add_axes([0.92, 0.02, 0.015, 0.43])
#    cbar = fig1.colorbar(image, cax=cax, label='dB')
#
#    plt.tight_layout()
#    plt.savefig(f'{out_dir}/{p}/{t_name}_{r_name}_{p}_beam_slices.png')
#    plt.close()


if __name__=='__main__':
    
    import argparse
    import numpy as np
    from pathlib import Path
    import concurrent.futures
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gs
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from scipy.stats import median_absolute_deviation as mad

    import sys
    sys.path.append('../decode_rf_data')
    from colormap import spectral
    
    # Custom spectral colormap
    cmap = spectral()
    
    parser = argparse.ArgumentParser(description="""
        Plot healpix map of reference data
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='../../outputs/tile_maps/compare_beams/',
            help='Output directory. Default=../../outputs/tile_maps/compare_beams/')
    parser.add_argument('--map_dir', metavar='\b', default='../../outputs/tile_maps/tile_maps_norm/',
            help='Tile map directory. Default=../../outputs/tile_maps/tile_maps_norm/')
    parser.add_argument('--fee_map', metavar='\b', default='../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz',
            help='Healpix FEE map of mwa tile. default=../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz')
    parser.add_argument('--nside', metavar='\b',type=int, default=32,help='Healpix Nside. Default = 32')
    
    
    args = parser.parse_args()
    
    out_dir     = Path(args.out_dir)
    map_dir     = Path(args.map_dir)
    fee_map     = args.fee_map
    nside       = args.nside
    
#    # make output dir if it doesn't exist
#    out_dir.mkdir(parents=True, exist_ok=True)
#
#    # MWA beam pointings
#    pointings = ['0', '2', '4', '41']
#
#    # Create output directory structure
#    for p in pointings:
#        Path(f'{out_dir}/{p}/').mkdir(parents=True, exist_ok=True)

    # find all map files
    map_files = [item for item in map_dir.glob('*.npz')]

    plt_grid(f'{map_dir}/S07XX_rf1XX_tile_maps.npz',f'{map_dir}/S07YY_rf1YY_tile_maps.npz' )

    # Parallization magic happens here
#    with concurrent.futures.ProcessPoolExecutor() as executor:
#        results = executor.map(beam_slice, map_files)
#    beam_slice(map_files[0])
