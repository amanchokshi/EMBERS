import sys
import unittest
import numpy as np
import skyfield as sf
from astropy.time import Time
from skyfield.api import Topos, load
sys.path.append('../../code/sat_ephemeris')
import sat_ephemeris as se


ts = load.timescale()
ts = load.timescale(builtin=True)

# Position of MWA site in Lat/Lon/Elevation
MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)


class Test_sat_ephemeris(unittest.TestCase):

    sats, epochs = se.load_tle('../../data/tests/TLE/25986.txt')
    epoch_range = se.epoch_ranges(epochs)
    t_arr, index_epoch = se.epoch_time_array(epoch_range, 0, 4)
    passes, alt, az = se.sat_pass(sats, t_arr, index_epoch)
    time_array, sat_alt, sat_az = se.ephem_data(t_arr, passes[0], alt, az)


    def test_load_tle(self):
        self.assertEqual(self.sats[0].model.satnum, 25986)
        self.assertEqual(self.epochs[0], 2458756.5)
        self.assertEqual(len(self.epochs), len(self.sats))


    def test_epoch_ranges(self):
        self.assertEqual(len(self.epoch_range), 20)
        self.assertEqual(self.epoch_range[0], 2458756.5)
    
    def test_epoch_time_array(self):
        self.assertEqual(self.t_arr[0].tt_calendar(), (2019, 9, 30, 0, 0, 0.0))
        self.assertEqual(self.index_epoch, 0)

    def test_sat_pass(self):
        self.assertEqual(self.passes[0][0], 1202)
        self.assertEqual(self.passes.shape, (3,2))
        self.assertEqual(self.alt.radians.shape, (10800,))
        self.assertEqual(self.alt.radians[0], -0.48374984876262395)
        self.assertEqual(self.az.radians.shape, (10800,))
        self.assertEqual(self.az.radians[0], 1.526615252001072)
        
    def test_ephem_data(self):
        self.assertEqual(len(self.sat_alt), len(self.sat_az))
        self.assertEqual(self.time_array[0], 1569806338.8159924)
        self.assertEqual(self.sat_alt[0], -1.150703355237408)
        self.assertEqual(self.sat_az[0], 3.8934109628585034)
    
    def test_sat_plot(self):
        plt = se.sat_plot(25986, self.sat_alt, self.sat_az, len(self.time_array), alpha=0.3)
        num_figures_before = plt.gcf().number
        plt = se.sat_plot(25986, self.sat_alt, self.sat_az, len(self.time_array), alpha=0.3)
        num_figures_after = plt.gcf().number
        self.assertEqual(num_figures_before+1, num_figures_after)

        

if __name__ == "__main__":
    unittest.main()
