import sys
import json
import unittest
from pathlib import Path
sys.path.append('../../code/beam_pointings')
import obs_pointings as op
import sort_pointings as sp


class Test_obs_pointings(unittest.TestCase):

    start_date  = '2019-10-01'
    stop_date   = '2019-10-10'
    time_zone   = 'Australia/Perth'
    out_dir     = './'
    meta_dir    = Path('../../data/tests/mwa_pointings')
    f_name      = Path('./ultimate_pointing_times.json')
    
    # create ultimate pointing.json
    start_gps, stop_gps, obs_length, pointings = sp.org_pointing_json(meta_dir)
    sp.combine_obs(start_gps, stop_gps, obs_length, pointings, out_dir) 
       
    # Test time_tree
    obs_time, obs_unix, obs_unix_end, obs_gps, obs_gps_end = op.time_tree(
            start_date,
            stop_date,
            time_zone)
        
    # Test read_json
    pointings, start_gps, stop_gps, obs_length = op.read_json(out_dir, f_name)
    
    
    def test_time_tree(self):
        self.assertEqual(self.obs_time[0], f'{self.start_date}-00:00')
        self.assertEqual(self.obs_gps[0], 1253894418.0)

    def test_read_json(self):
        self.assertEqual(self.pointings[0], 100)
        self.assertEqual(self.pointings[-1], 4)
        self.assertEqual(self.obs_length[0], 1056)

    def test_time_filter(self):
        occu = op.time_filter(self.start_gps[0], self.stop_gps[0], [self.obs_gps[0], self.obs_gps_end[0]])
        self.assertIsNone(occu)
        occu = op.time_filter(self.start_gps[12], self.stop_gps[12], [self.obs_gps[0], self.obs_gps_end[0]])
        self.assertEqual(occu,0.8255555555555556)

    def test_obs_pointing(self):
        
        op.obs_pointing(
                self.obs_time, 
                self.obs_gps, 
                self.obs_gps_end, 
                self.pointings,
                self.start_gps,
                self.stop_gps,
                self.start_date, 
                self.stop_date,
                self.out_dir
                )

        obs_json = Path('./obs_pointings.json')
        
        with open(obs_json, 'r') as oj:
            obs_point = json.load(oj)
            point_0 = obs_point["point_0"]

        self.assertTrue(obs_json.is_file())
        self.assertEqual(point_0[0], f'{self.start_date}-03:30')

        if obs_json.is_file():
            obs_json.unlink()
            self.f_name.unlink()


if __name__ == "__main__":
    unittest.main()


