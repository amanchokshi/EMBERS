import sys
import shutil
import unittest
import concurrent.futures
from itertools import repeat
from datetime import datetime, timedelta
sys.path.append('../../code/decode_rf_data')
import batch_waterfall as bf
from pathlib import Path
import rf_data as rf

# Test plt_single_waterfall
rf_file = 'S06XX_2019-10-01-14:30'
rf_dir = '../../data/rf_data/S06XX/2019-10-01'
bf.plt_single_waterfall(rf_dir, rf_file, './')

# Test waterfall_plot
tile_name = 'rf0XX'
out_dir = './'
data_dir = '../../data/rf_data'
date = '2019-10-01'
time_stamp = '2019-10-01-14:30'
bf.waterfall_plot(tile_name, out_dir, data_dir, date, time_stamp)



class Test_batch_waterfall(unittest.TestCase):

    
    def test_plt_single_waterfall(self):    
        self.assertTrue(Path(f'{rf_file}.png').is_file())

        if Path(f'{rf_file}.png').is_file() is True:
            Path(f'{rf_file}.png').unlink()

    def test_waterfall_plot(self):
        rf_path = Path(f'waterfalls/{date}/{time_stamp}/{tile_name}_{time_stamp}.png')
        self.assertTrue(rf_path.is_file())

        if rf_path.is_file() is True:
            shutil.rmtree(Path(f'waterfalls'))


if __name__ == "__main__":
    unittest.main()


