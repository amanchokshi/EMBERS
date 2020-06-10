import sys
import unittest

sys.path.append('../../code/decode_rf_data')
import rf_data as rf

power, times = rf.read_data('../../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
print(power.shape)
print(times.shape)

class Test_rf_data(unittest.TestCase):

    def test_read_data(self):
        power, times = rf.read_data('../../data/rf_data/S06XX/2019-10-01/S06XX_2019-10-01-14:30.txt')
        self.assertEqual(power.shape, (16655, 112))
        self.assertEqual(power.shape[0], times.shape[0])


