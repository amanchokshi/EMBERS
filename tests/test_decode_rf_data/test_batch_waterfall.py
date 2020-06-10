import sys
import unittest
import concurrent.futures
from datetime import datetime, timedelta
sys.path.append('../../code/decode_rf_data')
import batch_waterfall as bf
from pathlib import Path
import rf_data as rf

## Test plt_single_waterfall
#rf_file = 'S06XX_2019-10-01-14:30'
#rf_dir = '../../data/rf_data/S06XX/2019-10-01'
#bf.plt_single_waterfall(rf_dir, rf_file, './')

# Test waterfall_plot
data_dir = '../../data/rf_data'
start_date = '2019-10-01'
stop_date = '2019-10-01'
out_dir = './'

tiles = rf.tile_names()
t_start = datetime.strptime(start_date, '%Y-%m-%d')
t_stop = datetime.strptime(stop_date, '%Y-%m-%d')
n_days = (t_stop - t_start).days

dates = []
date_time = []

for i in range(n_days+1):
    day = t_start + timedelta(days=i)
    date = day.strftime('%Y-%m-%d')
    dates.append(date)
    d_t = []

    for j in range(48):
        t_delta = datetime.strptime(date,'%Y-%m-%d') + timedelta(minutes=30*j)
        d_time = t_delta.strftime('%Y-%m-%d-%H:%M')
        d_t.append(d_time)

        date_time.append(d_t)


data_dir = Path(data_dir)
out_dir = Path(out_dir)

for tile in tiles:
    for d in range(len(dates)):
        rf_path = data_dir/tile/dates[d]
        for t in date_time[d]:
            rf_name = f'{tile}_{t}'
            bf.waterfall_plot(t)
        #with concurrent.futures.ProcessPoolExecutor() as executor:
        #    results = executor.map(bf.waterfall_plot, date_time[d])

            #for result in results:
            #    print(result)

class Test_batch_waterfall(unittest.TestCase):

    
#    def test_plt_single_waterfall(self):    
#        self.assertTrue(Path(f'{rf_file}.png').is_file())
#
#        if Path(f'{rf_file}.png').is_file() is True:
#            Path(f'{rf_file}.png').unlink()

    def test_waterfall_plot(self):
        pass






if __name__ == "__main__":
    unittest.main()


