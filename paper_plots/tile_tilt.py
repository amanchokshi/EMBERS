import json

import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import spectral
#  from embers.tile_maps.beam_utils import plot_healpix
from matplotlib import pyplot as plt

_spec, _ = spectral()

nside = 256
data = "./interp_maps/pix_resi.json"

delta_az = 83.1912

p_0 = [0, 0]
p_2 = [90, 6.80880000000001]
p_4 = [270, 6.80880000000001]

with open(data) as f:
    tilt_resi = json.load(f)

# pixel of minimum residual power
phase_cen = {}

for tile in tilt_resi.keys():
    if "rf1" in tile and "_0" in tile:
        #  if "S06" in tile:
        min_pix = np.argmax(tilt_resi[tile])
        #  if min_pix <= 600:
        #  phase_cen[f"{tile}"] = np.degrees(np.array(hp.pix2ang(nside, min_pix)))
        phase_cen[f"{tile}"] = np.array(hp.pix2ang(nside, min_pix))

theta_XX = []
phi_XX = []
theta_YY = []
phi_YY = []
for i in list(phase_cen.keys()):
    if "XX" in i:
        theta_XX.append(phase_cen[i][0])
        phi_XX.append(phase_cen[i][1])
    else:
        theta_YY.append(phase_cen[i][0])
        phi_YY.append(phase_cen[i][1])

theta = np.mean([theta_XX, theta_YY], axis=0)
phi = np.mean([phi_XX, phi_YY], axis=0)


colors = _spec(np.linspace(0.17, 0.9, len(theta_XX)))

plt.style.use("seaborn")
figure = plt.figure(figsize=(6, 6))
ax = figure.add_subplot(111, polar=True)

ax.set_ylim(0, 1.8)
ax.set_rgrids([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6], angle=90)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

for i, c in enumerate(colors):
    plt.scatter(
        phi_YY[i],
        np.degrees(theta_YY[i]),
        marker="P",
        color=c,
        edgecolor="black",
        linewidth=0.7,
        s=77,
    )
    plt.scatter(
        phi_XX[i],
        np.degrees(theta_XX[i]),
        marker="X",
        color=c,
        edgecolor="black",
        linewidth=0.7,
        s=77,
    )
    plt.scatter(
        phi[i],
        np.degrees(theta[i]),
        marker="o",
        color=c,
        edgecolor="black",
        linewidth=0.7,
        s=77,
    )
#  plt.legend()
plt.tight_layout()
plt.show()
