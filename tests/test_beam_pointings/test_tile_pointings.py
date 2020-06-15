import sys
import json
import unittest
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append('../../code/sat_channels')
from window_chan_map import time_tree
sys.path.append('../../code/beam_pointings')
import tile_pointings as tip

class Test_sort_pointings(unittest.TestCase):

        
    with open(f'../../data/tests/obs_pointings.json') as table:
        data = json.load(table)
        start_date  = data['start_date']
        stop_date   = data['stop_date']
        point_0     = data['point_0']
        point_2     = data['point_2']
        point_4     = data['point_4']
        point_41    = data['point_41']

    data_dir = Path('../../data/rf_data')
        
    p_0, p_2, p_4, p_41 = tip.point_int('S06XX', start_date, stop_date, point_0, point_2, point_4, point_41, data_dir)
        
    def test_org_pointing_json(self):
        self.assertEqual(self.p_0, 0.5)
    
    def test_combine_obs(self):
        fig = plt.figure(figsize=(6,4))
        int_dict = {'S06XX':[7,6,5,4]}
        ax = tip.plt_pointing_hist(1, 1, 1, fig, int_dict=int_dict, tile='S06XX')
        plt.savefig('tile_point.png')

        point_png = Path('tile_point.png')
        self.assertTrue(point_png.is_file())

        if point_png.is_file():
            point_png.unlink()

if __name__ == "__main__":
    unittest.main()


