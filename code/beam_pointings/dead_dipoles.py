import numpy as np
from astropy.io import fits
from pathlib import Path

meta = Path('./1252145824.metafits')


hdu = fits.open(meta)
obsid = hdu[0].header['GPSTIME']
delays = hdu[1].data['Delays']

print(repr(hdu[0].header))

## Counts the number of flagged dipoles for each tile
#flags = []
#for i in range(delays.shape[0]):
#    flags.append(np.count_nonzero(delays[i] == 32))
#
#flags = np.array(flags)
## Seperates XX, YY tiles
#xx = flags[1::2]
#yy = flags[::2]
#xx_1 = np.count_nonzero(xx == 1)
#xx_2 = np.count_nonzero(xx == 2)
#yy_1 = np.count_nonzero(yy == 1)
#yy_2 = np.count_nonzero(yy == 2)
#metadata = [obsid, xx_1, xx_2, yy_1, yy_2]
#data.append(metadata)
#
tiles = [
        'HexS6',
        'HexS7',
        'HexS8',
        'HexS9',
        'HexS10',
        'HexS12',
        'HexS29',
        'HexS30',
        'HexS31',
        'HexS32',
        'HexS33',
        'HexS34',
        'HexS35',
        'HexS36'
        ]



