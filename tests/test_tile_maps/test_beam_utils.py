import shutil
from os import path
from pathlib import Path

import numpy as np
from embers.tile_maps.beam_utils import (chisq_fit_gain,
                                         healpix_cardinal_indices,
                                         healpix_cardinal_slices, map_slices,
                                         nan_mad, plot_healpix, plt_slice,
                                         poly_fit, rotate_map, chisq_fit_test)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))
clean_map = Path(f"{test_data}/tile_maps/tile_maps_clean/S06XX_rf0XX_tile_maps.npz")
map_data = np.load(clean_map, allow_pickle=True)
map_med = np.asarray([(np.nanmedian(i) if i != [] else np.nan) for i in map_data["0"]])

nside = 32


def test_rotate_map():
    map_rot = rotate_map(nside, angle=0, healpix_array=map_med)
    assert map_rot.all() == map_med.all()
