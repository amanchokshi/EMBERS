import sys
import unittest
import numpy as np
from pathlib import Path
sys.path.append('../../code/beam_pointings')
import sort_pointings as sp


class Test_sort_pointings(unittest.TestCase):

    meta_dir = Path('../../data/tests/mwa_pointings')
    out_dir = Path('./')

        
        
    def test_org_pointing_json(self):
        start_gps, stop_gps, obs_length, pointings = sp.org_pointing_json(self.meta_dir)
        
        self.assertEqual(len(start_gps),len(stop_gps))
        self.assertEqual(start_gps[1], stop_gps[0])
        self.assertEqual(np.median(obs_length), 120)
        self.assertEqual(pointings[0], 100) 
        self.assertEqual(pointings[-1], 4) 

    
    def test_combine_obs(self):
        start_gps, stop_gps, obs_length, pointings = sp.org_pointing_json(self.meta_dir)
        sp.combine_obs(start_gps, stop_gps, obs_length, pointings, self.out_dir) 
        
        ulti_json = Path('ultimate_pointing_times.json')
        self.assertTrue(ulti_json.is_file())

        if ulti_json.is_file():
            ulti_json.unlink()

if __name__ == "__main__":
    unittest.main()


