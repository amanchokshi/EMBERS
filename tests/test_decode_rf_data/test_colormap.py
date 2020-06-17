import sys
import unittest

sys.path.append("../../code/decode_rf_data")
import colormap as cmap


class Test_colormap(unittest.TestCase):

    # Test first and last colors of custom colormaps
    def test_spec(self):
        self.assertEqual(cmap.spectral().colors[0][0], 0.14453125)
        self.assertEqual(cmap.spectral().colors[-1][0], 0.8555759803921569)

    # Test first and last colors of custom colormaps
    def test_jade(self):
        self.assertEqual(cmap.jade()[0].colors[0][0], 0.08048525330056805)
        self.assertEqual(cmap.jade()[0].colors[-1][0], 0.8736808458053672)

    # Test first and last colors of custom colormaps
    def test_kelp(self):
        self.assertEqual(cmap.kelp()[0].colors[0][0], 0.06813277130228468)
        self.assertEqual(cmap.kelp()[0].colors[-1][0], 0.8376025159364093)


if __name__ == "__main__":
    unittest.main()
