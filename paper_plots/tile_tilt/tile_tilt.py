import json

import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import jade, spectral
#  from embers.tile_maps.beam_utils import plot_healpix
from matplotlib import pyplot as plt

#  from tilt_fit import reproject_map

_spec, _ = spectral()
_jade, _ = jade()

nside = 256
data = "./interp_maps/pix_resi.json"
out_dir = "./interp_maps/beam_offsets"
m_256 = "./interp_maps/0"
m_32 = "../embers_out/tile_maps/tile_maps/tile_maps_clean"

ref_p = 243.17
ref_t = 0.476


with open(data) as f:
    tilt_resi = json.load(f)

# pixel of minimum residual power
phase_cen = {}

for tile in tilt_resi.keys():
    if "rf1" in tile and "_0" in tile:
        min_pix = np.argmax(tilt_resi[tile])
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

#  ax.set_ylim(0, 1.5)
#  ax.set_rgrids([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], angle=90)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

for i, c in enumerate(colors):
    #  plt.scatter(
    #  (phi_YY[i]),
    #  np.degrees(theta_YY[i]),
    #  marker="P",
    #  color=c,
    #  edgecolor="black",
    #  linewidth=0.7,
    #  s=77,
    #  )
    #  plt.scatter(
    #  (phi_XX[i]),
    #  np.degrees(theta_XX[i]),
    #  marker="X",
    #  color=c,
    #  edgecolor="black",
    #  linewidth=0.7,
    #  s=77,
    #  )
    plt.scatter(
        (phi[i]),
        np.degrees(theta[i]),
        marker="s",
        color=c,
        edgecolor="black",
        linewidth=0.7,
        s=84,
        alpha=0.9,
    )
#  plt.legend()
ax.set_xticklabels(["N", "", "E", "", "S", "", "W", ""])
plt.tight_layout()
plt.savefig("./interp_maps/beam_offsets/polar_offsets.png")
plt.close()


#  m_06_32 = "../../embers_out/tile_maps/tile_maps/tile_maps_clean/S08XX_rf1XX_tile_maps.npz"
#  m_06_32 = np.load(m_06_32, allow_pickle=True)["0"]
#  m_06_32 = np.asarray([(np.nanmedian(j) if j != [] else np.nan) for j in m_06_32])

#  m_06_256 = "./interp_maps/0/S08XX_rf0XX_0_N256.npz"
#  m_06_256 = np.load(m_06_256)["healpix"]

#  fee_256 = "./interp_maps/FEE/fee_XX_0_N256.npz"
#  fee_256 = np.load(fee_256)["healpix"]

#  res_256 = m_06_256 - fee_256
#  res_256[np.where(fee_256 < -30)] = np.nan

#  pix = np.argmax(tilt_resi["S08XX_rf1XX_0"])
#  rot_256 = reproject_map(256, pix, healpix_array=m_06_256)
#  res_256_rot = rot_256 - fee_256
#  res_256_rot[np.where(fee_256 < -30)] = np.nan


#  fig1 = plt.figure(figsize=(8.0, 8.0))
#  fig1.suptitle(f"Pix: {pix} (θ, ɸ): {np.degrees(hp.pix2ang(256, pix))}")

#  plt.subplot(2, 2, 1)
#  plot_healpix(data_map=m_06_32, vmin=-50, vmax=0, cmap=_jade, cbar=False)
#  plt.subplot(2, 2, 2)
#  plot_healpix(data_map=m_06_256, vmin=-50, vmax=0, cmap=_jade, cbar=False)
#  plt.subplot(2, 2, 3)
#  plot_healpix(data_map=res_256, vmin=-4, vmax=4, cmap="RdYlGn", cbar=False)
#  plt.subplot(2, 2, 4)
#  plot_healpix(data_map=res_256_rot, vmin=-4, vmax=4, cmap="RdYlGn", cbar=False)

#  plt.tight_layout()
#  plt.savefig(f"{out_dir}/S08XX_rf1_0.png")
