import shutil
from os import path
from pathlib import Path

from embers.mwa_utils.mwa_dipoles import (download_metafits, find_flags,
                                          mwa_flagged_dipoles)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


def test_download_metafits():
    download_metafits(1, 0, f"{test_data}/mwa_utils")
    meta = Path(f"{test_data}/mwa_utils/mwa_metafits/1253960760.metafits")
    assert meta.is_file()
    if meta.is_file():
        shutil.rmtree(f"{test_data}/mwa_utils/mwa_metafits")


def test_find_flags():
    download_metafits(2, 0, f"{test_data}/mwa_utils")
    find_flags(f"{test_data}/mwa_utils")
    plt = Path(f"{test_data}/mwa_utils/flagged_dipoles.png")
    assert plt.is_file()
    if plt.is_file():
        plt.unlink()
    meta = Path(f"{test_data}/mwa_utils/mwa_metafits/1253960760.metafits")
    assert meta.is_file()
    if meta.is_file():
        shutil.rmtree(f"{test_data}/mwa_utils/mwa_metafits")


def test_mwa_flagged_dipoles():
    mwa_flagged_dipoles(1, f"{test_data}/mwa_utils", 0)
    plt = Path(f"{test_data}/mwa_utils/flagged_dipoles.png")
    assert plt.is_file()
    if plt.is_file():
        plt.unlink()
    meta = Path(f"{test_data}/mwa_utils/mwa_metafits/1253960760.metafits")
    assert meta.is_file()
    if meta.is_file():
        shutil.rmtree(f"{test_data}/mwa_utils/mwa_metafits")
