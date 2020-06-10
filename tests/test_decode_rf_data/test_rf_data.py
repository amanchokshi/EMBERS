import sys
import unittest

sys.path.append('../../code/decode_rf_data')
import rf_data as rf

#power, times = rf.read_data('../../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
#print(power.shape)
#print(times.shape)

class Test_rf_data(unittest.TestCase):

    def test_read_data(self):
        power, times = rf.read_data('../../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
        
        self.assertEqual(power.shape[0], times.shape[0])
        self.assertEqual(power.shape, (16655, 112))

    def test_tile_names(self):
            tiles = [
                    'rf0XX', 'rf0YY', 'rf1XX', 'rf1YY',
                    'S06XX', 'S06YY', 'S07XX', 'S07YY',
                    'S08XX', 'S08YY', 'S09XX', 'S09YY',
                    'S10XX', 'S10YY', 'S12XX', 'S12YY',
                    'S29XX', 'S29YY', 'S30XX', 'S30YY',
                    'S31XX', 'S31YY', 'S32XX', 'S32YY',
                    'S33XX', 'S33YY', 'S34XX', 'S34YY',
                    'S35XX', 'S35YY', 'S36XX', 'S36YY'
                    ]
            self.assertEqual(rf.tile_names(), tiles)


if __name__ == "__main__":
    unittest.main()


