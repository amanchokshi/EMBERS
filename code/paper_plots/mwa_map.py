import json
import argparse
import numpy as np
from pathlib import Path
from astropy.io import fits
import matplotlib.pylab as pl
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="""
        Look for dead dipoles
        """)

parser.add_argument(
        '--metafits', metavar='\b', default='./../../outputs/beam_pointings/metafits/1262518040.metafits', 
        help='Directory where metafits files live. Default=./../../outputs/beam_pointings/metafits/1262518040.metafits')

parser.add_argument(
        '--out_dir', metavar='\b', default='./../../outputs/beam_pointings/', 
        help='Directory where metafits files live. Default=./../../outputs/paper_plots/')

args = parser.parse_args()
metafits     = Path(args.metafits)
out_dir       = Path(args.out_dir)

# list of all tiles used in this experiment
tiles_14 = [
        'HexS6', 'HexS7',
        'HexS8', 'HexS9',
        'HexS10', 'HexS12',
        'HexS29', 'HexS30',
        'HexS31', 'HexS32',
        'HexS33', 'HexS34', 
        'HexS35', 'HexS36'
        ]



tiles = []
N = []
E = []

hdu = fits.open(metafits)
tile_names = hdu[1].data['TileName']
north = hdu[1].data['North']
east = hdu[1].data['East']
pols = hdu[1].data['Pol']

for t in range(len(tile_names)):
    #if 'HexS' in tile_names[t] and 'X' in pols[t]: 
    if 'X' in pols[t]: 
        tiles.append(tile_names[t])
        N.append(north[t])
        E.append(east[t])

rf_n = [-39.74, -10.88]
rf_e = [36.26, 67.48]

n_14 = []
e_14 = []
for n in tiles_14:
    idx_n = np.where(np.array(tiles) == n)[0][0]
    e_14.append(E[idx_n])
    n_14.append(N[idx_n])
    
#print(n_14)
#print(e_14)


nice_fonts = {
        # Use LaTeX to write all text
        "text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 10,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        }

plt.rcParams.update(nice_fonts)


fig, ax = plt.subplots(figsize=(3.6, 3.6))
#fig, ax = plt.subplots(figsize=(16, 16))
ax.scatter(E, N, marker='.', color='#f9ed69')
ax.scatter(rf_e, rf_n, marker='.', color='red')
ax.scatter(e_14, n_14, marker='.', color='blue')

#for i, txt in enumerate(tiles):
#    ax.annotate(txt, (E[i], N[i]))

ax.set_aspect('equal')
ax.set_ylim(-100, 400)
ax.set_xlim(-200, 300)
#ax.set_ylim(50, 200)
#ax.set_xlim(-50, 100)
ax.set_ylabel('North [m]')
ax.set_xlabel('East [m]')

plt.savefig('test.png', transparent=True, bbox_inches='tight')

