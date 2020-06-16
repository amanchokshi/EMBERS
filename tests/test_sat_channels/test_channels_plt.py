import sys
import json
import shutil
import unittest
import numpy as np
from pathlib import Path
sys.path.append('../../code/sat_channels')
import channels_plt as cp
import window_chan_map as wcm
sys.path.append('../../code/decode_rf_data')
from colormap import spectral
cmap = spectral()


class Test_channels_plt(unittest.TestCase):
    
    # Window chan argparse variables
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
    
    # Test time filter
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

    
    def test_plt_waterfall_pass(self):
        cp.plt_waterfall_pass(
                self.power, 25414, self.intvl[0], 
                self.intvl[1], self.start_date, 
                cmap, chs=[11,41,46,52],
                good_ch=11, out_dir='./')

        waterfall = Path('2019-10-01_25414_[11, 41, 46, 52]_waterfall.png')
        self.assertTrue(waterfall.is_file())
        if waterfall.is_file:
            waterfall.unlink()


    def test_plt_channel_basic(self):
        cp.plt_channel_basic(
                './', self.times, self.power[:, 11],
                11, min(self.power[:, 11]), 
                max(self.power[:, 11]), 
                self.noise_t, self.pow_thresh, 
                25414, self.start_date)
        ch_plt = Path('2019-10-01_25414_11_channel.png')
        self.assertTrue(ch_plt.is_file())
        if ch_plt.is_file():
            ch_plt.unlink()

    def test_sat_plot(self):
        cp.sat_plot('./', ['A', 'B', 'C'], 25414, [[1], [2], [3]], [[1],[2],[3]], 3, self.start_date, 'polar')
        polar = Path('2019-10-01_25414_polar.png')
        self.assertTrue(polar.is_file())
        if polar.is_file():
            polar.unlink()


if __name__ == "__main__":
    unittest.main()
