import itertools
from pathlib import Path

import healpy as hp
import numpy as np
from embers.rf_tools.colormaps import jade
from healpy.rotator import euler_matrix_new as euler
from matplotlib import pyplot as plt
from scipy import interpolate

jade, _ = jade()

fee_map = "../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz"
nside = 32


#  def reproject_map(nside, pixel=None, healpix_array=None):
def reproject_map(nside, healpix_array=None):
    """Reproject beam map to be centered around healpix pixel instead of zenith.

    :param nside: Healpix nside
    :param pixel: Healpix index of new phase centre
    :param healpix_array: Input healpix array to be rotated

    :returns:
        - Reprojected healpix map

    """

    # coordinated of new phase centre
    #  theta, phi = hp.pix2ang(nside, pixel)
    theta, phi = [1, 0]

    vec = hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))
    eu_mat = euler(-phi, -theta, phi, deg=True)
    rot_map = hp.rotator.rotateVector(eu_mat, vec)
    new_hp_inds = hp.vec2pix(nside, rot_map[0], rot_map[1], rot_map[2])

    return healpix_array[new_hp_inds]


def plot_healpix(
    data_map=None,
    rot=(0, 90, 180),
    fig=None,
    sub=None,
    title=None,
    vmin=None,
    vmax=None,
    cmap=None,
    cbar=True,
):
    """Yeesh do some healpix magic to plot the thing

    :param data_map: Healpix input map to plot
    :param fig: Figure number to use
    :param sub: Matplotlib subplot syntax
    :param title: Plot title
    :param vmin: Colormap minimum
    :param vmax: Colormap maximum
    :param cmap: Matplotlib :class:`~matplotlib.colors.ListedColormap`
    :param cbar: If True, plot a colorbar

    :returns:
        - Plot of healpix map

    """

    # Disable cryptic healpy warnings. Can't figure out where they originate
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    hp.delgraticules()
    hp.orthview(
        map=data_map,
        coord="E",
        fig=fig,
        half_sky=True,
        rot=rot,
        xsize=1200,
        title=title,
        sub=sub,
        min=vmin,
        max=vmax,
        cmap=cmap,
        notext=True,
        hold=True,
        cbar=cbar,
        return_projected_map=True,
    )

    hp.graticule(dpar=10, coord="E", color="k", alpha=0.7, dmer=45, lw=0.4, ls=":")

    # Altitude grid
    hp.projtext(
        00.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "0",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )
    hp.projtext(
        30.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "30",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )
    hp.projtext(
        60.0 * (np.pi / 180.0),
        225.0 * (np.pi / 180),
        "60",
        color="k",
        coord="E",
        fontsize=6,
        fontweight="light",
    )

    # NSEW
    hp.projtext(
        80.0 * (np.pi / 180.0),
        000.0 * (np.pi / 180.0),
        r"$N  $",
        coord="E",
        color="w",
        fontweight="light",
        verticalalignment="top",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        090.0 * (np.pi / 180.0),
        r"$E  $",
        coord="E",
        color="w",
        fontweight="light",
        horizontalalignment="right",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        180.0 * (np.pi / 180.0),
        r"$S  $",
        coord="E",
        color="w",
        fontweight="light",
        verticalalignment="bottom",
    )
    hp.projtext(
        80.0 * (np.pi / 180.0),
        270.0 * (np.pi / 180.0),
        r"$W  $",
        coord="E",
        color="w",
        fontweight="light",
        horizontalalignment="left",
    )


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


#  maps = beam_maps(
#  "../embers_out/tile_maps/tile_maps/tile_maps_clean/S06XX_rf0XX_tile_maps.npz"
#  )[0][0]

#  fee_m = np.load(fee_map, allow_pickle=True)
#  fee = fee_m["0"][0]

#  plot_healpix(data_map=maps, vmin=-50, vmax=0)
#  plt.savefig("rot_i.png")
#  plt.close()

#  residuals = maps - fee
#  residuals[np.where(fee < -30)] = np.nan
#  residuals[np.where(maps == np.nan)] = np.nan

#  plot_healpix(data_map=residuals, vmin=-5, vmax=5)
#  plt.savefig("rot_i_res.png")
#  plt.close()

#  #  map_rot = reproject_map(nside, 0, healpix_array=maps)
#  map_rot = reproject_map(nside, healpix_array=maps)

#  residuals = map_rot - fee
#  residuals[np.where(fee < -30)] = np.nan
#  residuals[np.where(map_rot == np.nan)] = np.nan

#  plot_healpix(data_map=map_rot, vmin=-50, vmax=0)
#  plt.savefig("rot_ii.png")
#  plt.close()

#  plot_healpix(data_map=residuals, vmin=-5, vmax=5)
#  plt.savefig("rot_ii_res.png")
#  plt.close()


def polyfit2d(x, y, z, order=3):
    ncols = (order + 1) ** 2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order + 1), range(order + 1))
    for k, (i, j) in enumerate(ij):
        G[:, k] = x ** i * y ** j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m


def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order + 1), range(order + 1))
    z = np.zeros_like(x)
    for a, (i, j) in zip(m, ij):
        z += a * x ** i * y ** j
    return z


z_map = beam_maps(
    "../embers_out/tile_maps/tile_maps/tile_maps_clean/S06XX_rf0XX_tile_maps.npz"
)[0][0]

pix_ind = np.arange(hp.nside2npix(nside))
x, y, _ = hp.pix2vec(nside, pix_ind)
z = z_map[:6144]
mask = np.isnan(z)
x = x[:6144][~mask]
y = y[:6144][~mask]
z = z[~mask]


def circular_mask(h, w):

    center = (int(w/2), int(h/2))
    radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask


xnew_edges, ynew_edges = np.mgrid[-1:1:1024j, -1:1:1024j]
xnew = xnew_edges[:-1, :-1] + np.diff(xnew_edges[:2, 0])[0] / 2.0
ynew = ynew_edges[:-1, :-1] + np.diff(ynew_edges[0, :2])[0] / 2.0

tck = interpolate.bisplrep(x, y, z)
znew = interpolate.bisplev(xnew[:, 0], ynew[0, :], tck)
mask = circular_mask(znew.shape[0], znew.shape[1])
znew[~mask] = np.nan
plt.figure()
plt.pcolormesh(xnew_edges, ynew_edges, znew, shading="flat", vmin=-50, vmax=0)
plt.colorbar()
plt.title("Interpolated function.")
plt.show()

#  # Fit a 3rd order, 2d polynomial
#  m = polyfit2d(x, y, z, order=5)

#  # Evaluate it on a grid...
#  nx, ny = 1024, 1024
#  xx, yy = np.meshgrid(
#  np.linspace(x.min(), x.max(), nx), np.linspace(y.min(), y.max(), ny)
#  )
#  zz = polyval2d(xx, yy, m)

#  # Plot
#  plt.imshow(zz, extent=(x.min(), y.max(), x.max(), y.min()))
#  #  plt.scatter(x, y, c=z, s=1)
#  plt.show()

#  #  plt.scatter(x[:6144], y[:6144], c=z, vmin=-50, vmax=0)
#  #  ax = plt.gca()
#  #  ax.set_aspect(1)
#  #  plt.colorbar()
#  #  plt.show()
