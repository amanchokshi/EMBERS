import sys
import shutil
import argparse
import unittest
import numpy as np
from pathlib import Path

sys.path.append("../../code/align_data")
sys.path.append("../../code/decode_rf_data")
from batch_align import save_aligned


class Test_batch_align(unittest.TestCase):

    ref = "rf0XX"
    aut = "S06XX"
    ref_dir = "../../data/rf_data/rf0XX/2019-10-01"
    aut_dir = "../../data/rf_data/S06XX/2019-10-01"
    savgol_window_1 = 11
    savgol_window_2 = 15
    polyorder = 2
    interp_type = "cubic"
    interp_freq = 1
    out_dir = "./"
    date = "2019-10-01"
    timestamp = "2019-10-01-14:30"

    def test_save_aligned(self):
        save_aligned(
            self.ref,
            self.aut,
            self.ref_dir,
            self.aut_dir,
            self.savgol_window_1,
            self.savgol_window_2,
            self.polyorder,
            self.interp_type,
            self.interp_freq,
            self.out_dir,
            self.date,
            self.timestamp,
        )

        aligned_npz = Path(
            f"{self.date}/{self.timestamp}/{self.ref}_{self.aut}_{self.timestamp}_aligned.npz"
        )
        aligned_data = np.load(aligned_npz)
        self.assertTrue(aligned_npz.is_file())
        self.assertTrue(
            aligned_data["ref_p_aligned"].shape, aligned_data["tile_p_aligned"].shape
        )

        if aligned_npz.is_file() is True:
            shutil.rmtree(Path(f"{self.date}"))


if __name__ == "__main__":
    unittest.main()
