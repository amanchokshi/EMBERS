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

for t in tiles:
    idx = np.where(tile_names == t)
    print(tile_names[idx])



