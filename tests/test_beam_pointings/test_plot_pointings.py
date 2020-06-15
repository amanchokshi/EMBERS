import sys
import unittest
import numpy as np
from pathlib import Path
sys.path.append('../../code/beam_pointings')
import sort_pointings as sp
import plot_pointings as pp


class Test_sort_pointings(unittest.TestCase):

    meta_dir = Path('../../data/tests/mwa_pointings')
    out_dir = Path('./')
    f_name = 'ultimate_pointing_times.json'
    int_thresh = 2

    start_gps, stop_gps, obs_length, pointings = sp.org_pointing_json(meta_dir)
    sp.combine_obs(start_gps, stop_gps, obs_length, pointings, out_dir)
        
    time_point, point = pp.point_integration(f_name, out_dir, int_thresh)
        
    def test_point_integration(self):
        self.assertEqual(self.time_point, [24.273333333333333, 3.4444444444444446, 6.995555555555556])
        self.assertEqual(self.point, [0, 33, 135])

        if Path(self.f_name).is_file:
            Path(self.f_name).unlink()

    def test_pointing_hist(self):
        pp.pointing_hist(self.time_point, self.point, self.out_dir)

        self.assertTrue(Path('pointing_integration.png').is_file())

        if Path('pointing_integration.png').is_file():
            Path('pointing_integration.png').unlink()

    #def test_combine_obs(self):
    #    
    #    ulti_json = Path('ultimate_pointing_times.json')
    #    self.assertTrue(ulti_json.is_file())

    #    if ulti_json.is_file():
    #        ulti_json.unlink()

if __name__ == "__main__":
    unittest.main()


