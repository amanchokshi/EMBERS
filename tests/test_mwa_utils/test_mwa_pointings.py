import shutil
from os import path
from pathlib import Path

from embers.mwa_utils.mwa_pointings import (clean_meta_json, combine_pointings,
                                            download_meta, mwa_point_meta,
                                            obs_pointings, plt_hist_array,
                                            point_integration, pointing_hist,
                                            rf_obs_times, tile_integration)

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


def test_download_meta():
    download_meta(
        "2019-10-01", "2019-10-10", 1, f"{test_data}/mwa_utils/mwa_meta_tmp", wait=0
    )
    meta = Path(f"{test_data}/mwa_utils/mwa_meta_tmp/mwa_pointings/page_001.json")
    assert meta.is_file()
    if meta.is_file():
        shutil.rmtree(f"{test_data}/mwa_utils/mwa_meta_tmp")

