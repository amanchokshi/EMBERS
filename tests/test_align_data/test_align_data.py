import sys
import unittest

sys.path.append("../../code/align_data")
sys.path.append("../../code/decode_rf_data")
from align_data import savgol_interp
from rf_data import read_data


class Test_align_data(unittest.TestCase):
    def test_savgol_interp(self):

        (
            ref_t,
            ref_p,
            tile_t,
            tile_p,
            ref_p_aligned,
            tile_p_aligned,
            time_array,
        ) = savgol_interp(
            "../../data/rf_data/rf0XX/2019-10-10/rf0XX_2019-10-10-02:30.txt",
            "../../data/rf_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt",
            savgol_window_1=11,
            savgol_window_2=15,
            polyorder=2,
            interp_type="cubic",
            interp_freq=1,
        )

        self.assertEqual(ref_p_aligned.shape, tile_p_aligned.shape)
        self.assertEqual(time_array.shape[0], 1778)


if __name__ == "__main__":
    unittest.main()
