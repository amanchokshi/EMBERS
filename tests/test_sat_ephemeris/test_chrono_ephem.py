import sys
import json
import unittest
from pathlib import Path
sys.path.append('../../code/sat_ephemeris')
import chrono_ephem as ce


            
class Test_chrono_ephem(unittest.TestCase):

    def test_obs_times(self):
        start_date = '2019-10-10'
        stop_date = '2019-10-10'
        time_zone = 'Australia/Perth'

        obs_times, obs_unix, obs_unix_end = ce.obs_times(time_zone, start_date, stop_date)
        self.assertEqual(obs_times[0], '2019-10-10-00:00')
        self.assertEqual(len(obs_unix), len(obs_unix_end))
        self.assertEqual(obs_unix_end[0] - obs_unix[0], 1800)

    def test_write_json(self):
        l = [1,2,3,4]
        data = {'test':l}
        ce.write_json(data, filename='test.json', out_dir='./')
        
        with open('test.json', 'r') as data:
            data_r = json.load(data)
            data_json = data_r['test']
        
        self.assertEqual(data_json, l)
        Path('test.json').unlink()
    
    def test_interp_ephem(self):
        
        with open('../../data/tests/ephem_json/25986.json') as ephem:
            sat_ephem = json.load(ephem)
            
            # Extract data from json dictionary
            t_array = sat_ephem['time_array']
            s_alt   = sat_ephem['sat_alt']
            s_az    = sat_ephem['sat_az']
        
            # here, we're looping over each satellite pass with a single sat ephem file 
            # to check which observation window it falls in
            
            # Don't consider passes with less than 4 samples (~ 1 min = 3*20s)
            time_interp, sat_alt, sat_az = ce.interp_ephem(
                    t_array[0],
                    s_alt[0],
                    s_az[0],
                    'cubic',
                    1)

        self.assertEqual(len(sat_alt), len(sat_az))
        self.assertEqual(len(time_interp), len(sat_az))
        self.assertEqual(time_interp[0], 1569806331.0)
        self.assertEqual(sat_alt[0], -1.5292196979547075)
        self.assertEqual(sat_az[0], 3.9006136945178773)


if __name__ == "__main__":
    unittest.main()
