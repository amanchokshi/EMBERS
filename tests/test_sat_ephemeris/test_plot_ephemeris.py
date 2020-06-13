import sys
import json
import unittest
sys.path.append('../../code/sat_ephemeris')
import plot_ephemeris as pe


class Test_plot_ephemeris(unittest.TestCase):
    
    json_dir = './'
        
    def test_plot_ephem(self):    

        pe.plot_ephem(self.tle_dir, self.cadence, self.out_dir, self.sat_id)
        out_file = Path(f'./{self.sat_id}.json')

        self.assertTrue(out_file.is_file())
        self.assertEqual(sat_id_json[0], self.sat_id)
        out_file.unlink()


if __name__ == "__main__":
    unittest.main()
