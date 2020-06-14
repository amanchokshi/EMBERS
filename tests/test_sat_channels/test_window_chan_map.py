import sys
import json
import unittest
import numpy as np
sys.path.append('../../code/sat_channels')
import window_chan_map as wcm
sys.path.append('../../code/decode_rf_data')
from colormap import spectral


class Test_window_chan_map(unittest.TestCase):

    start_date = '2019-10-01'
    stop_date = '2019-10-01'
    ali_npz = '../../data/tests/align_data/2019-10-01/2019-10-01-14:30/rf0XX_S06XX_2019-10-01-14:30_aligned.npz'
    chrono_file = '../../data/tests/chrono_json/2019-10-01-14:30.json'
    sat_thresh = 1
    noise_thresh = 3
    pow_thresh = 15
    occ_thresh = 0.8

    dates, date_time = wcm.time_tree(start_date, stop_date)
    power, times = wcm.read_aligned(ali_file=ali_npz)
    power, noise_t = wcm.noise_floor(sat_thresh, noise_thresh, data=power)
    
    with open(chrono_file) as chrono:
        chrono_ephem = json.load(chrono)
        norad_list = []
        for s in range(len(chrono_ephem)):
            norad_list.append(chrono_ephem[s]["sat_id"][0])
        
        norad_index = norad_list.index(norad_list[0])
        norad_ephem = chrono_ephem[norad_index]
            
        rise_ephem  = norad_ephem["time_array"][0] 
        set_ephem   = norad_ephem["time_array"][-1]
        sat_alt     = norad_ephem["sat_alt"]

        intvl = wcm.time_filter(rise_ephem, set_ephem, np.asarray(times))

    good_chan = wcm.good_chans(
            ali_npz, chrono_file, 
            norad_list[0], sat_thresh, 
            noise_thresh, pow_thresh, 
            occ_thresh, '2019-10-01', 
            '2019-10-01-14:30', False)

    def test_time_tree(self):
        self.assertEqual(self.start_date, self.dates[0])
        self.assertEqual(self.date_time[0][0], f'{self.start_date}-00:00')
        self.assertEqual(self.date_time[0][-1], f'{self.start_date}-23:30')
        self.assertEqual(len(self.date_time[0]), 48)


    def test_read_aligned(self):    
        self.assertEqual(self.power.shape[0], self.times.shape[0])
        self.assertEqual(self.power.shape[0], 1779)
        self.assertEqual(self.power.shape[-1], 112)


    def test_noise_floor(self):
        self.assertEqual(self.noise_t, 0.8962653625488282)
        self.assertTrue(np.isclose(0, np.median(self.power), atol=1e-3))


    def test_time_filter(self):
        self.assertEqual(self.intvl, [27, 621])
#
#    def test_good_chans(self):
#        pass
#
#    def test_window_chan_map(self):
#        pass

        

if __name__ == "__main__":
    unittest.main()
