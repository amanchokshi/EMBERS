from pathlib import Path

import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import jade
from embers.tile_maps.beam_utils import plot_healpix
from matplotlib import pyplot as plt
from scipy import interpolate
from scipy import interpolate as interp
from scipy.spatial.transform import Rotation as R

jade, _ = jade()

fee_map = "../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz"
nside = 32


def beam_maps(f):
    """Returns pairs of [AUT,FEE] maps, for each pointing"""
    f_name, _ = Path(f).name.split(".")
    t_name, r_name, _, _ = f_name.split("_")

    # pointings = ['0','2','4','41']
    pointings = ["0", "2", "4"]

    maps = []

    # load data from map .npz file
    tile_map = np.load(f, allow_pickle=True)
    fee_m = np.load(fee_map, allow_pickle=True)

    for p in pointings:

        tile = tile_map[p]

        if "XX" in t_name:
            fee = fee_m[p][0]
        else:
            fee = fee_m[p][1]

        # Visualize the tile map and diff map
        # healpix meadian map
        tile_med = np.asarray([(np.nanmedian(j) if j != [] else np.nan) for j in tile])

        residuals = tile_med - fee
        residuals[np.where(fee < -30)] = np.nan
        residuals[np.where(tile_med == np.nan)] = np.nan

        maps.append([tile_med, residuals])

    return maps


b_map = beam_maps(
    "../embers_out/tile_maps/tile_maps/tile_maps_clean/S06XX_rf0XX_tile_maps.npz"
)[0][0]

fee_m = np.load(fee_map, allow_pickle=True)
fee = fee_m["0"][0]


def beam_spline(beam_map):

    # Input data
    pix_ind = np.arange(hp.nside2npix(nside))
    x, y, _ = hp.pix2vec(nside, pix_ind)
    z = beam_map[:6144]
    mask = np.isnan(z)
    x = x[:6144][~mask]
    y = y[:6144][~mask]
    z = z[~mask]

    #  def circular_mask(h, w):

    #  center = (int(w / 2), int(h / 2))
    #  radius = min(center[0], center[1], w - center[0], h - center[1])

    #  Y, X = np.ogrid[:h, :w]
    #  dist_from_center = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2)

    #  mask = dist_from_center <= radius
    #  return mask

    #  xnew_edges, ynew_edges = np.mgrid[-1:1:512j, -1:1:512j]
    #  xnew = xnew_edges[:-1, :-1] + np.diff(xnew_edges[:2, 0])[0] / 2.0
    #  ynew = ynew_edges[:-1, :-1] + np.diff(ynew_edges[0, :2])[0] / 2.0

    tck = interpolate.bisplrep(x, y, z, s=5000)
    #  znew = interpolate.bisplev(xnew[:, 0], ynew[0, :], tck)
    #  mask = circular_mask(znew.shape[0], znew.shape[1])
    #  znew[~mask] = np.nan

    #  return (xnew, ynew, znew)
    znew = interpolate.bisplev([1, 0, -1], [0, 1, 0], tck)
    print(znew)


#  fee_x, fee_y, fee_z = beam_spline(fee)
#  beam_x, beam_y, beam_z = beam_spline(b_map)

#  res_z = beam_z - fee_z
#  fee_z[np.isnan(fee_z)] = -80
#  mask = np.where(fee_z < -30)
#  res_z[mask] = np.nan

#  # stack the meshgrid to position vectors
#  xyz = np.vstack([fee_x.ravel(), fee_y.ravel(), fee_z.ravel(),]).T

#  r = R.from_euler("zyz", [45, 10, -10], degrees=True)

#  rot_beam = r.apply(xyz)
#  rot_beam = rot_beam.T

#  x_r = rot_beam[0].reshape(511, 511)
#  y_r = rot_beam[0].reshape(511, 511)
#  z_r = rot_beam[0].reshape(511, 511)

#  plt.imshow(z_r)
#  plt.show()
#  plt.figure()
#  #  plt.pcolormesh(beam_x, beam_y, res_z, shading="flat", vmin=-5, vmax=5)
#  plt.pcolormesh(x_r, y_r, fee_z, shading="flat")
#  plt.colorbar()
#  plt.title("Interpolated function.")
#  plt.show()


def beam_rbf(beam_map):

    # Input data
    i_pix = np.arange(hp.nside2npix(nside))
    theta, phi = hp.pix2ang(nside, i_pix)
    power = beam_map[:6144]
    mask = np.isnan(power)
    theta = theta[:6144][~mask]
    phi = phi[:6144][~mask]
    power = power[~mask]
    rbfi = interp.Rbf(theta, phi, power, function="cubic")

    t, p = hp.pix2ang(128, np.arange(int(hp.nside2npix(128))))

    p_128 = rbfi(t, p)
    #  plot_healpix(data_map=p_128, vmin=-50, vmax=0)
    #  plt.savefig("rbf.png")
    return p_128


b_128 = beam_rbf(b_map)
f_128 = beam_rbf(fee)

res_128 = b_128 - f_128
mask = np.where(f_128 < -30)
res_128[mask] = np.nan

plot_healpix(data_map=res_128, vmin=-5, vmax=5)
plt.savefig("res_rbf.png")
