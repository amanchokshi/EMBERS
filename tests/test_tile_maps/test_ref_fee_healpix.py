import shutil
from os import path
from pathlib import Path

import pkg_resources
from embers.tile_maps.ref_fee_healpix import create_model, ref_healpix_save

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


def test_create_model_first():
    refXX = pkg_resources.resource_filename(
        "embers.kindle", "data/ref_models/MWA_reference_tile_FarField_XX.ffe"
    )
    healpix_XX, theta_mesh, phi_mesh, power_XX, theta = create_model(
        32, file_name=refXX
    )
    assert round(healpix_XX[0]) == 0


def test_create_model_last():
    refXX = pkg_resources.resource_filename(
        "embers.kindle", "data/ref_models/MWA_reference_tile_FarField_XX.ffe"
    )
    healpix_XX, theta_mesh, phi_mesh, power_XX, theta = create_model(
        32, file_name=refXX
    )
    assert round(healpix_XX[-1]) == -2567420


def test_ref_healpix_save():
    ref_healpix_save(32, f"{test_data}/tile_maps/tmp")
    npz = Path(f"{test_data}/tile_maps/tmp/ref_dipole_models.npz")
    assert npz.is_file()
    if npz.is_file():
        shutil.rmtree(f"{test_data}/tile_maps/tmp")
