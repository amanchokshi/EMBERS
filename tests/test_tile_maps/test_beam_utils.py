from os import path
from pathlib import Path

import numpy as np
from embers.tile_maps.beam_utils import (chisq_fit_gain, chisq_fit_test,
                                         healpix_cardinal_indices,
                                         healpix_cardinal_slices, map_slices,
                                         nan_mad, plt_slice,
                                         poly_fit, rotate_map)
from matplotlib import pyplot as plt

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))
clean_map = Path(f"{test_data}/tile_maps/tile_maps_clean/S06XX_rf0XX_tile_maps.npz")
map_data = np.load(clean_map, allow_pickle=True)
map_med = np.asarray([(np.nanmedian(i) if i != [] else np.nan) for i in map_data["0"]])

nside = 32


def test_rotate_map():
    map_rot = rotate_map(nside, angle=np.pi, healpix_array=map_med)
    assert map_rot[1] == map_med[3]


def test_rotate_map_flip():
    map_rot = rotate_map(nside, angle=0, healpix_array=map_med, flip=True)
    assert map_rot[1] == map_med[0]


def test_rotate_map_save():
    npz = Path(f"{test_data}/tile_maps/tmp_rot_map.npz")
    rotate_map(nside, angle=0, healpix_array=map_med, flip=True, savetag=npz)
    assert npz.is_file()
    if npz.is_file():
        npz.unlink()


def test_healpix_cardinal_indices():
    NS, EW = healpix_cardinal_indices(nside)
    assert [NS[0], EW[0]] == [5968, 6000]


def test_healpix_cardinal_slices():
    NS, EW = healpix_cardinal_slices(nside, map_med, 90)
    assert [round(np.nanmedian(NS[0])), round(np.nanmedian(EW[0]))] == [-37.0, -24.0]


def test_nan_mad():
    map_mad = nan_mad(map_data["0"])
    assert map_mad[0] == 0


def test_map_slices():
    NS, EW = map_slices(nside, map_data["0"], 90)
    assert round(NS[2][0]) == -89


def test_poly_fit():
    fit = poly_fit(
        np.linspace(1, 10, 10), np.linspace(3, 13, 10), np.linspace(1, 10, 10), 3
    )
    assert round(fit[0]) == 3


def test_chisq_fit_test():
    pval = chisq_fit_test(
        data=np.linspace(1, 10, 10), model=np.linspace(3, 13, 10), offset=20
    )
    assert round(pval) == 1


def test_chisq_fit_gain():
    gain = chisq_fit_gain(np.linspace(1, 10, 10), np.linspace(3, 13, 10))
    assert round(gain[0]) == -3


def test_plt_slice():
    NS, _ = map_slices(nside, map_data["0"], 90)
    fig = plt.figure()
    ax = plt_slice(
        fig=fig,
        sub=[1, 1, 1],
        zen_angle=NS[2],
        map_slice=NS[0],
        model_slice=NS[0],
        delta_pow=NS[0],
        slice_label="test",
        pow_fit=NS[0],
        xlabel=True,
        ylabel=True,
    )
    assert type(ax).__name__ == "AxesSubplot"


def test_plt_slice_labels():
    NS, _ = map_slices(nside, map_data["0"], 90)
    fig = plt.figure()
    ax = plt_slice(
        fig=fig,
        sub=[1, 1, 1],
        zen_angle=NS[2],
        map_slice=NS[0],
        model_slice=NS[0],
        delta_pow=NS[0],
        slice_label="test",
        pow_fit=NS[0],
        xlabel=False,
        ylabel=False,
    )
    assert type(ax).__name__ == "AxesSubplot"
