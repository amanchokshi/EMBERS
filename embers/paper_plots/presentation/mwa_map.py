import json
import argparse
import numpy as np
from pathlib import Path
from astropy.io import fits
import matplotlib.pylab as pl
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(
    description="""
        Look for dead dipoles
        """
)

parser.add_argument(
    "--metafits",
    metavar="\b",
    default="./../../../outputs/beam_pointings/metafits/1262518040.metafits",
    help="Directory where metafits files live. Default=./../../../outputs/beam_pointings/metafits/1262518040.metafits",
)

parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./../../../outputs/paper_plots/",
    help="Directory where metafits files live. Default=./../../../outputs/paper_plots/",
)

args = parser.parse_args()
metafits = Path(args.metafits)
out_dir = Path(args.out_dir)

# list of all tiles used in this experiment
tiles_14 = [
    "HexS6",
    "HexS7",
    "HexS8",
    "HexS9",
    "HexS10",
    "HexS12",
    "HexS29",
    "HexS30",
    "HexS31",
    "HexS32",
    "HexS33",
    "HexS34",
    "HexS35",
    "HexS36",
]


# Read metafits file and extract positions of each tile in N, E coordinates
tiles = []
N = []
E = []
hdu = fits.open(metafits)
tile_names = hdu[1].data["TileName"]
north = hdu[1].data["North"]
east = hdu[1].data["East"]
pols = hdu[1].data["Pol"]

for t in range(len(tile_names)):
    if "X" in pols[t]:
        tiles.append(tile_names[t])
        N.append(north[t])
        E.append(east[t])

# Positions of the references
rf_n = [-39.74, -10.88]
rf_e = [36.26, 67.48]

# positions of 14 tiles used
n_14 = []
e_14 = []
for n in tiles_14:
    idx_n = np.where(np.array(tiles) == n)[0][0]
    e_14.append(E[idx_n])
    n_14.append(N[idx_n])


nice_fonts = {
    # Use LaTeX to write all text
    # "text.usetex": True,
    "font.family": "sans-serif",
    # Use 10pt font in plots, to match 10pt font in document
    "axes.labelsize": 10,
    "font.size": 10,
    # Make the legend/label fonts a little smaller
    "legend.fontsize": 6,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "ytick.color": "#696969",
    "xtick.color": "#696969",
    "axes.labelcolor": "#696969",
    "axes.edgecolor": "#696969",
}

plt.rcParams.update(nice_fonts)


fig, ax = plt.subplots(figsize=(3.6, 3.6))
# fig, ax = plt.subplots(figsize=(16, 16))
ax.scatter(E, N, marker="s", s=10.0, color="orange", label="MWA")
ax.scatter(e_14, n_14, marker="s", s=10.0, color="cornflowerblue", label="AUT")
ax.scatter(rf_e, rf_n, marker="s", s=10.0, color="crimson", label="REF")
# ax.scatter(-8.5, 0.9, marker='s', s=10.0, color='green', label='Center')

for i, txt in enumerate([r"$Ref_1$", r"$Ref_0$"]):
    ax.annotate(
        txt,
        (rf_e[i], rf_n[i]),
        xytext=(1.28 * rf_e[i], rf_n[i]),
        va="center",
        color="w",
        fontsize=7,
    )

leg = ax.legend(
    loc="lower right",
    frameon=True,
    framealpha=0.5,
    markerscale=1,
    prop={"weight": "heavy"},
)
leg.get_frame().set_edgecolor("w")
leg.get_frame().set_facecolor("#222222")
for text in leg.get_texts():
    text.set_color("white")
for l in leg.legendHandles:
    l.set_alpha(1)

# ax.legend()
ax.set_aspect("equal")
ax.set_ylim(-100, 400)
ax.set_xlim(-200, 300)
ax.set_ylabel("North [m]")
ax.set_xlabel("East [m]")

plt.savefig(f"mwa_map.pdf", transparent=True, bbox_inches="tight")
