import sys
import unittest
sys.path.append('../../code/beam_pointings')
import date_gps as gps


class Test_date_gps(unittest.TestCase):
        
    def test_convert_gps(self):    
        self.assertEqual(gps.convert_gps('2019-10-10'), 1254700818)

if __name__ == "__main__":
    unittest.main()


