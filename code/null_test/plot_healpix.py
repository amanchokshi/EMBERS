# Healpix plotting script adapted from Dr. Jack Line's code
# https://github.com/JLBLine/MWA_ORBCOMM

import numpy as np
import healpy as hp
from scipy.stats import median_absolute_deviation as mad


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
    ref, _, _ = f_name.split('_')

    # Plot beam map
    map_data = np.load(f, allow_pickle=True)
    ref_map = map_data['ref_map']
    ref_counter = map_data['ref_counter']
    
    # compute the median for every pixel array
    ref_map_med = [(np.mean(i) if i != [] else np.nan ) for i in ref_map]
    vmin = np.nanmin(ref_map_med)
    vmax = np.nanmax(ref_map_med)

    fig = plt.figure(figsize=(8,10))
    fig.suptitle(f'Reference Beam Healpix: {ref}', fontsize=16)
    plot_healpix(data_map=np.asarray(ref_map_med),sub=(1,1,1), cmap=cmap, vmin=vmin, vmax=vmax)

    plt.savefig(f'{out_dir}/{f_name}.png',bbox_inches='tight')
    plt.close()


    # Plot MAD
    # compute the median for every pixel array
    ref_map_med = [(np.median(i) if i != [] else np.nan ) for i in ref_map]

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


    vmin = np.nanmin(ref_map_mad)
    vmax = np.nanmax(ref_map_mad)

    fig = plt.figure(figsize=(8,10))
    fig.suptitle(f'Healpix MAD: {ref}', fontsize=16)
    plot_healpix(data_map=np.asarray(ref_map_mad),sub=(1,1,1), cmap=cmap, vmin=vmin, vmax=vmax)

    plt.savefig(f'{out_dir}/{f_name}_mad.png',bbox_inches='tight')
    plt.close()

    
    # Plot counts in pix
    fig = plt.figure(figsize=(8,10))
    fig.suptitle(f'Healpix Pixel Counts: {ref}', fontsize=16)
    plot_healpix(data_map=np.asarray(ref_counter),sub=(1,1,1), cmap=cmap, vmin=0, vmax=2400)

    plt.savefig(f'{out_dir}/{f_name}_counts.png',bbox_inches='tight')
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
    from colormap import spectral
    
    # Custom spectral colormap
    cmap = spectral()
    

    parser = argparse.ArgumentParser(description="""
        Plot healpix map of reference data
        """)
    
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/null_test/',help='Output directory. Default=./../../outputs/null_test/')
    
    args = parser.parse_args()
    
    out_dir = Path(args.out_dir)

    map_files = [item for item in out_dir.glob('*.npz')]
    print(map_files)

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(map_plots, map_files)


