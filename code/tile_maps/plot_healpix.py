# Healpix plotting script adapted from Dr. Jack Line's code
# https://github.com/JLBLine/MWA_ORBCOMM

import sys
import numpy as np
import healpy as hp
from scipy.stats import median_absolute_deviation as mad

sys.path.append('../sat_ephemeris')
from sat_ids import norad_ids

def plot_healpix(data_map=None,sub=None,title=None,vmin=None,vmax=None,cmap=None):
    '''Yeesh do some healpix magic to plot the thing'''
    
    # Disable cryptic healpy warnings. Can't figure out where they originate
    import warnings
    warnings.filterwarnings("ignore", category=RuntimeWarning) 
    
    if vmin == None:
        if cmap == None:
            half_sky = hp.orthview(
                    map=data_map,coord='E',
                    half_sky=True,xsize=400,
                    title=title,rot=(0,90,0),
                    sub=sub,notext=True, return_projected_map=True)
        else:
            half_sky = hp.orthview(
                    map=data_map,coord='E',
                    half_sky=True,xsize=400,
                    title=title,rot=(0,90,0), sub=sub,cmap=cmap,
                    notext=True,return_projected_map=True)
    else:
        if cmap == None:
            half_sky = hp.orthview(
                    map=data_map,coord='E'
                    ,half_sky=True,xsize=400,rot=(0,90,0),
                    title=title,sub=sub,min=vmin,max=vmax,
                    notext=True,return_projected_map=True)
        else:
            half_sky = hp.orthview(
                    map=data_map,coord='E',
                    half_sky=True,xsize=400,rot=(0,90,0),
                    title=title,sub=sub,min=vmin,max=vmax,
                    cmap=cmap,notext=True,return_projected_map=True)

    hp.graticule(dpar=10,coord='E',color='k',alpha=0.3,dmer=45)
   
    # Altitude grid
    hp.projtext(0.0*(np.pi/180.0), 0.0, '0', coord='E')
    hp.projtext(30.0*(np.pi/180.0), 0.0, '30', coord='E')
    hp.projtext(60.0*(np.pi/180.0), 0.0, '60', coord='E')

    # Azimuth grid
    hp.projtext(90.0*(np.pi/180.0), 00.0*(np.pi/180.0), r'$0^\circ$', coord='E',color='k',verticalalignment='top', fontsize=12)
    hp.projtext(90.0*(np.pi/180.0), 90.0*(np.pi/180.0), r'$90^\circ$', coord='E',color='k',horizontalalignment='right', fontsize=12)
    hp.projtext(90.0*(np.pi/180.0), 180.0*(np.pi/180.0), r'$180^\circ$', coord='E',color='k', fontsize=12)
    hp.projtext(90.0*(np.pi/180.0), 270.0*(np.pi/180.0), r'$270^\circ$', coord='E',color='k', fontsize=12)
    
    # NSEW 
    hp.projtext(90.0*(np.pi/180.0), 045.0*(np.pi/180.0), r'$N  $', coord='E',color='k',verticalalignment='top', horizontalalignment='right', fontsize=14)
    hp.projtext(90.0*(np.pi/180.0), 135.0*(np.pi/180.0), r'$E  $', coord='E',color='k',horizontalalignment='right', fontsize=14)
    hp.projtext(90.0*(np.pi/180.0), 225.0*(np.pi/180.0), r'$S  $', coord='E',color='k', fontsize=14)
    hp.projtext(90.0*(np.pi/180.0), 315.0*(np.pi/180.0), r'$W  $', coord='E',color='k', verticalalignment='top', horizontalalignment='left', fontsize=14)


def map_plots(f):
    
    f_name, _ = f.name.split('.')
    tile, ref, _, _ = f_name.split('_')
   
    pointings = ['0','2','4']

    # load data from map .npz file
    map_data = np.load(f, allow_pickle=True)
    map_data = {key:map_data[key].item() for key in map_data}
    tile_maps = map_data['healpix_maps']
    #tile_counters = map_data['healpix_counters']

    for p in pointings:

        tile_map = tile_maps[p]
        #tile_counter = tile_counters[p]
        tile_counter = [len(m) for m in tile_map]
    
        # Plot BEAM
        # compute the median for every pixel array
        tile_map_mean = [(np.mean(i) if i != [] else np.nan ) for i in tile_map]
        tile_map_mean_scaled = np.asarray([(i - np.nanmax(tile_map_mean[:5000])) for i in tile_map_mean])
        vmin = np.nanmin(tile_map_mean_scaled)
        #vmax = np.nanmax(tile_map_mean_scaled)
        vmax = 0

        fig = plt.figure(figsize=(8,10))
        fig.suptitle(f'Healpix Map: {tile}/{ref} @ {p}', fontsize=16)
        plot_healpix(data_map=np.asarray(tile_map_mean_scaled),sub=(1,1,1), cmap=jade, vmin=vmin, vmax=vmax)
        plt.savefig(f'{out_dir}/tile_maps/{tile}_{ref}_{p}_map.png',bbox_inches='tight')
        plt.close()
           
        # Plot MAD 
        tile_map_mad = []
        for j in tile_map:
            if j != []:
                j = np.asarray(j)
                j = j[~np.isnan(j)]
                tile_map_mad.append(mad(j))
            else:
                tile_map_mad.append(np.nan)

        tile_map_mad = np.asarray(tile_map_mad)
        
        tile_map_mad[np.where(tile_map_mad == np.nan)] = np.nanmean(tile_map_mad)


        vmin = np.nanmin(tile_map_mad)
        vmax = np.nanmax(tile_map_mad)

        fig = plt.figure(figsize=(8,10))
        fig.suptitle(f'Healpix MAD: {tile}/{ref} @ {p}', fontsize=16)
        plot_healpix(data_map=np.asarray(tile_map_mad),sub=(1,1,1), cmap=jade, vmin=vmin, vmax=vmax)
        plt.savefig(f'{out_dir}/tile_errors/{tile}_{ref}_{p}_mad.png',bbox_inches='tight')
        plt.close()


        # Plot counts in pix
            
        fig = plt.figure(figsize=(8,10))
        fig.suptitle(f'Healpix Pixel Counts: {tile}/{ref} @ {p}', fontsize=16)
        plot_healpix(data_map=np.asarray(tile_counter),sub=(1,1,1), cmap=jade, vmin=0, vmax=400)

        plt.savefig(f'{out_dir}/tile_counts/{tile}_{ref}_{p}_counts.png',bbox_inches='tight')
        plt.close()


def good_maps(f):
    
    f_name, _ = f.name.split('.')
    tile, ref, _, _ = f_name.split('_')
   
    pointings = ['0','2','4']

    # load data from map .npz file
    tile_data = np.load(f, allow_pickle=True)
    tile_data = {key:map_data[key].item() for key in tile_data}
    ratio_map = tile_data['ratio_map']  
    
    # Good sats from which to make plots
    good_sats = [
            25338, 25984, 25985,
            28654, 40086, 40087,
            40091, 41179, 41180,
            41182, 41183, 41184,
            41185, 41187, 41188,
            41189, 44387
            ]
    
    for p in pointings:
        
        # Empty good map
        good_map = [[] for pixel in range(hp.nside2npix(nside))]
        
        # append to good map from all good sat data
        for sat in good_sats:
            for pix in range(hp.nside2npix(nside)):
                good_map[pix].extend(ratio_map[p][sat][pix])
        
        good_map_med = [(np.median(i) if i != [] else np.nan ) for i in good_map]
        good_map_scaled = np.asarray([(i - np.nanmax(good_map_med[:5000])) for i in good_map_med])

        
        fig = plt.figure(figsize=(8,10))
        fig.suptitle(f'Good Map: {tile}/{ref} @ {p}', fontsize=16)
        plot_healpix(data_map=good_map_scaled, sub=(1,1,1), cmap=jade, vmin=np.nanmin(good_map_scaled), vmax=0)
        plt.savefig(f'{out_dir}/good_maps/{p}/tile_maps/{tile}_{ref}_{p}_good_map.png',bbox_inches='tight')
        plt.close()


        # Plot MAD 
        good_map_mad = []
        for j in good_map:
            if j != []:
                j = np.asarray(j)
                j = j[~np.isnan(j)]
                good_map_mad.append(mad(j))
            else:
                good_map_mad.append(np.nan)

        good_map_mad = np.asarray(tile_map_mad)
        
        good_map_mad[np.where(good_map_mad == np.nan)] = np.nanmedian(good_map_mad)

        vmin = np.nanmin(tile_map_mad)
        vmax = np.nanmax(tile_map_mad)

        fig = plt.figure(figsize=(8,10))
        fig.suptitle(f'Good Map MAD: {tile}/{ref} @ {p}', fontsize=16)
        plot_healpix(data_map=np.asarray(tile_map_mad),sub=(1,1,1), cmap=jade, vmin=vmin, vmax=vmax)
        plt.savefig(f'{out_dir}/good_maps/{p}/tile_errors/{tile}_{ref}_{p}_good_map_errors.png',bbox_inches='tight')
        plt.close()


        # Plot counts in pix
            
        good_map_counts = [len(i) for i in good_map]
        
        fig = plt.figure(figsize=(8,10))
        fig.suptitle(f'Good Map Counts: {tile}/{ref} @ {p}', fontsize=16)
        plot_healpix(data_map=np.asarray(tile_counter),sub=(1,1,1), cmap=jade, vmin=0, vmax=300)
        plt.savefig(f'{out_dir}/good_maps/{p}/tile_counts/{tile}_{ref}_{p}_good_map_counts.png',bbox_inches='tight')
        plt.close()


def sat_maps(f)
    
    f_name, _ = f.name.split('.')
    tile, ref, _, _ = f_name.split('_')
   
    pointings = ['0','2','4']

    # list of all possible satellites
    sat_ids = list(norad_ids.values())
    
    # load data from map .npz file
    tile_data = np.load(f, allow_pickle=True)
    tile_data = {key:map_data[key].item() for key in tile_data}
    ratio_map = tile_data['ratio_map']  
    
    for p in pointings:
        for s in sat_ids:

            fig = plt.figure(figsize=(8,10))
            fig.suptitle(f'Satellite [{s}]: {tile}/{ref} @ {p}', fontsize=16)
            ratio_sat_med = [(np.median(i) if i != [] else np.nan ) for i in tile_data['ratio_map'][p][s]]
            ratio_sat_scaled = np.asarray([(i - np.nanmax(ratio_sat_med[:5000])) for i in ratio_sat_med])
            plot_healpix(data_map=ratio_sat_scaled, sub=(1,1,1), cmap=jade)
            plt.savefig(f'{out_dir}/sat_maps/{p}/{s}_{p}_{tile}_{ref}_passes.png',bbox_inches='tight')
            plt.close()



if __name__=='__main__':
    
    import argparse
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from pathlib import Path
    import concurrent.futures
    
    import sys
    sys.path.append('../decode_rf_data')
    from colormap import spectral, jade, kelp
    
    # Custom spectral colormap
    spec = spectral()
    jade, _ = jade()
    kelp, _ = kelp()
    

    parser = argparse.ArgumentParser(description="""
        Plot healpix map of reference data
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/tile_maps/',help='Output directory. Default=./../../outputs/tile_maps/')
    parser.add_argument('--map_dir', metavar='\b', default='./../../outputs/tile_maps/tile_maps_sats',help='Output directory. Default=./../../outputs/tile_maps/tile_maps_sats')
    
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)
    map_dir = Path(args.map_dir)
    
    pointings = ['0','2','4']
    map_files = [item for item in map_dir.glob('*.npz')]

    # Create output directory tree
    for p in pointings:
        Path(f'{out_dir}/good_maps/{p}/tile_maps/').mkdir(parents=True, exist_ok=True)
        Path(f'{out_dir}/good_maps/{p}/tile_counts/').mkdir(parents=True, exist_ok=True)
        Path(f'{out_dir}/good_maps/{p}/tile_errors/').mkdir(parents=True, exist_ok=True)
        Path(f'{out_dir}/sat_maps/{p}/').mkdir(parents=True, exist_ok=True)

     # Depreciated: plotting maps from all data   
#    # Parallization magic happens here
#    with concurrent.futures.ProcessPoolExecutor() as executor:
#        results = executor.map(map_plots, map_files)

#    # Parallization magic happens here
#    with concurrent.futures.ProcessPoolExecutor() as executor:
#        results = executor.map(good_maps, map_files)

    # plot maps for all sats, for only one tile_ref pair
    # S08XX has good time coverage
    sat_maps(Path(f'{map_dir}/S08XX_rf0XX_sat_maps.npz'))

