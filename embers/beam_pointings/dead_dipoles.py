import json
import numpy as np
from pathlib import Path
from astropy.io import fits
import matplotlib.pylab as pl
import matplotlib.pyplot as plt


def find_flag(metafits, tiles):

    # dictionary to save flag info.
    # 0 indicates no dead dipole
    # 1-16 dipole index
    flags = {}
    flags['obsid'] = []
    for t in tiles:
        for p in ['X', 'Y']:
            flags[f'{t}{p}'] = []

    meta_files = [item for item in metafits.glob('*.metafits')]
    for m in meta_files: 

        hdu = fits.open(m)
        mode = hdu[0].header['MODE']
        obsid = hdu[0].header['GPSTIME']
        delays = hdu[1].data['Delays']
               
        tile_names = hdu[1].data['TileName']
        delays = hdu[1].data['Delays']
        pols = hdu[1].data['Pol']
        
        if mode == 'HW_LFILES':

            flags['obsid'].append(int(obsid))
            
            for t in tiles:
                idx = np.where(tile_names == t)
                for i in range(2):
                    idx_p = idx[0][i]
                    
                    t_name = tile_names[idx_p]
                    t_pol = pols[idx_p]
                    t_del = delays[idx_p]
                    t_flag = np.where(t_del == 32)[0] + 1
                     
                    if t_flag.size != 0:
                        flags[f'{t_name}{t_pol}'].append(int(t_flag[0]))
                    else:
                        flags[f'{t_name}{t_pol}'].append(0)

        hdu.close()

    return flags



if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="""
            Look for dead dipoles
            """)
    
    parser.add_argument(
            '--metafits', metavar='\b', default='./../../outputs/beam_pointings/metafits/', 
            help='Directory where metafits files live. Default=./../../outputs/beam_pointings/metafits/')
    
    parser.add_argument(
            '--out_dir', metavar='\b', default='./../../outputs/beam_pointings/', 
            help='Directory where metafits files live. Default=./../../outputs/beam_pointings/')
    
    args = parser.parse_args()
    metafits     = Path(args.metafits)
    out_dir       = Path(args.out_dir)
    
    
    # list of all tiles used in this experiment
    tiles = [
            'HexS6', 'HexS7',
            'HexS8', 'HexS9',
            'HexS10', 'HexS12',
            'HexS29', 'HexS30',
            'HexS31', 'HexS32',
            'HexS33', 'HexS34', 
            'HexS35', 'HexS36']

    flags = find_flag(metafits, tiles)

    #for m in meta_files:
    #    find_flag(m)
    
    # Save flagging info to json file
    with open(f'{out_dir}/flags.json', 'w') as outfile:
        json.dump(flags, outfile, indent=4)
    
    keys = list(flags.keys())
    
    n = len(keys) - 1
    colors = pl.cm.Spectral(np.linspace(0,1,n))
    
    plt.style.use('seaborn')
    
    fig, axs = plt.subplots(4,7, figsize=(18, 9), sharex=True, sharey=True,)
    axs = axs.ravel()
    
    for i in range(n):
        axs[i].scatter(flags['obsid'], flags[keys[i+1]], color=colors[i], linewidths=0.1, s=28 ,alpha=0.88, label=keys[i+1])
        axs[i].set_ylim(-1, 17)
        axs[i].set_yticks([0, 4, 8, 12, 16])
        
        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.2)
    
        # place a text box in upper left in axes coords
        axs[i].text(0.67, 0.93, f'{keys[i+1]}', transform=axs[i].transAxes, fontsize=10,
            verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    plt.savefig(f'{out_dir}/dead_dipoles.png')
    
    
    
