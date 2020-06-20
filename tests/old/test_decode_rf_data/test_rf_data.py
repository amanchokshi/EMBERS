import sys
import unittest

sys.path.append("../../code/decode_rf_data")
import rf_data as rf


class Test_rf_data(unittest.TestCase):

    power, times = rf.read_data(
        "../../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt"
    )

    # Test the shape of data and that the time & power arrays have the same length
    def test_read_data(self):
        self.assertEqual(self.power.shape[0], self.times.shape[0])
        self.assertEqual(self.power.shape, (16655, 112))

    # Test tile names. Hardcoded in rf_data
    def test_tile_names(self):
        tiles = [
            "rf0XX",
            "rf0YY",
            "rf1XX",
            "rf1YY",
            "S06XX",
            "S06YY",
            "S07XX",
            "S07YY",
            "S08XX",
            "S08YY",
            "S09XX",
            "S09YY",
            "S10XX",
            "S10YY",
            "S12XX",
            "S12YY",
            "S29XX",
            "S29YY",
            "S30XX",
            "S30YY",
            "S31XX",
            "S31YY",
            "S32XX",
            "S32YY",
            "S33XX",
            "S33YY",
            "S34XX",
            "S34YY",
            "S35XX",
            "S35YY",
            "S36XX",
            "S36YY",
        ]
        self.assertEqual(rf.tile_names(), tiles)

    # Test that a matplotlib plot figure is created
    def test_plot_waterfall(self):
        import matplotlib.pyplot as plt

        num_figures_before = plt.gcf().number
        plt = rf.plot_waterfall(self.power, self.times, None)
        num_figures_after = plt.gcf().number
        self.assertEqual(num_figures_before + 1, num_figures_after)


if __name__ == "__main__":
    unittest.main()
