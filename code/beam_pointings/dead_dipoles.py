import numpy as np
from astropy.io import fits
from pathlib import Path

meta = Path('./1252145824.metafits')


hdu = fits.open(meta)
obsid = hdu[0].header['GPSTIME']
delays = hdu[1].data['Delays']

tile_names = hdu[1].data['TileName']
delays = hdu[1].data['Delays']
pols = hdu[1].data['Pol']

tiles = [
        'HexS6', 'HexS7',
        'HexS8', 'HexS9',
        'HexS10', 'HexS12',
        'HexS29', 'HexS30',
        'HexS31', 'HexS32',
        'HexS33', 'HexS34', 
        'HexS35', 'HexS36'
        ]

flags = {}

for t in tiles:
    for p in ['X', 'Y']:
        flags[f'{t}{p}'] = []


for t in tiles:
    idx = np.where(tile_names == t)
    for i in range(2):
        idx_p = idx[0][i]
        
        t_name = tile_names[idx_p]
        t_pol = pols[idx_p]
        t_del = delays[idx_p]
        t_flag = np.where(t_del == 32)[0] + 1
         
        if t_flag.size != 0:
            print(t_name, t_pol, t_del, t_flag[0])
            flags[f'{t_name}{t_pol}'].append(t_flag[0])
        else:
            print(t_name, t_pol, t_del, 0)
            flags[f'{t_name}{t_pol}'].append(0)

print(flags)




