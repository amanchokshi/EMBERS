import sys
import unittest
sys.path.append('../../code/sat_ephemeris')
import sat_ids


class Test_sat_ids(unittest.TestCase):

    def test_norad_ids(self):    

        norad_dict = sat_ids.norad_ids

        self.assertEqual(len(norad_dict.keys()), 73)
        
        

if __name__ == "__main__":
    unittest.main()


