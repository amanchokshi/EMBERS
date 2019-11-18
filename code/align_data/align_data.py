import sys
sys.path.append('../decode_rf_data')

import rf_data as rf


r_power, r_times = rf.read_data('./../../data/rf0XX_2019-10-10-02:30.txt')

tile_power, tile_times = rf.read_data('./../../data/S10XX_2019-10-10-02:30.txt')

print((r_times[1]))
print((tile_times[1]))

