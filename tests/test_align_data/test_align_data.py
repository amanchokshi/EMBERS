import sys
import unittest
sys.path.append('../../code/align_data')
from align_data import savgol_interp


class Test_align_data(unittest.TestCase):

    def test_savgol_interp(self):    
        self.assertEqual(cmap.spectral().colors[0][0], 0.14453125)
        self.assertEqual(cmap.spectral().colors[-1][0], 0.8555759803921569)
        

if __name__ == "__main__":
    unittest.main()


