import json
import argparse
import numpy as np
from pathlib import Path
from astropy.io import fits

parser = argparse.ArgumentParser(description="""
        Look for dead dipoles
        """)

parser.add_argument(
        '--metafits', metavar='\b', default='./../../outputs/beam_pointings/metafits/', 
        help='Directory where metafits files live. Default=./../../outputs/beam_pointings/metafits/')

parser.add_argument(
        '--output', metavar='\b', default='./../../outputs/beam_pointings/flags.json', 
        help='Directory where metafits files live. Default=./../../outputs/beam_pointings/flags.json')

args = parser.parse_args()
metafits     = Path(args.metafits)
output       = Path(args.output)

meta_files = [item for item in metafits.glob('*.metafits')]

# list of all tiles used in this experiment
tiles = [
        'HexS6', 'HexS7',
        'HexS8', 'HexS9',
        'HexS10', 'HexS12',
        'HexS29', 'HexS30',
        'HexS31', 'HexS32',
        'HexS33', 'HexS34', 
        'HexS35', 'HexS36'
        ]

# dictionary to save flag info.
# 0 indicates no dead dipole
# 1-16 dipole index
flags = {}
for t in tiles:
    for p in ['X', 'Y']:
        flags[f'{t}{p}'] = []


def find_flag(meta):
    hdu = fits.open(meta)
    obsid = hdu[0].header['GPSTIME']
    delays = hdu[1].data['Delays']
    
    tile_names = hdu[1].data['TileName']
    delays = hdu[1].data['Delays']
    pols = hdu[1].data['Pol']
    
    
    
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

for m in meta_files:
    find_flag(m)

print(flags)


with open(f'{output}', 'w') as outfile:
    json.dump(flags, outfile, indent=4)


