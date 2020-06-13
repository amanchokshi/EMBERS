import sys
import json
import unittest
from pathlib import Path
from skyfield.api import Topos, load
sys.path.append('../../code/sat_ephemeris')
import batch_sat_ephem as ba



class Test_batch_sat_ephem(unittest.TestCase):

    # Position of MWA site in Lat/Lon/Elevation
    MWA = Topos(latitude=-26.703319, longitude=116.670815, elevation_m=337.83)
    
    tle_dir = './../../data/tests/TLE/'
    cadence = 360
    out_dir = './'
    sat_id = 25986

    def test_sat_json(self):    
        ba.sat_json(self.tle_dir, self.cadence, self.out_dir, self.sat_id)
        out_file = Path(f'./{self.sat_id}.json')

        with open(out_file, 'r') as ephem:
            sat_ephem = json.load(ephem)
            sat_id_json = sat_ephem['sat_id']
            altitude = sat_ephem['sat_alt']


        self.assertTrue(out_file.is_file())
        self.assertEqual(sat_id_json[0], self.sat_id)
        self.assertEqual(altitude[0][0],-7.013352046042244)
        out_file.unlink()
        

if __name__ == "__main__":
    unittest.main()

