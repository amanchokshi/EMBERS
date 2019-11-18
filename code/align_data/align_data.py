import sys
sys.path.append('../decode_rf_data')

import rf_data as rf

tile_list = rf.tile_names()


for i in rf.tile_names():
    _, times = rf.read_data(f'/media/achokshi/satellites/tiles_data/{i}/2019-11-15/{i}_2019-11-15-03:30.txt')
    print(f'{i}: {len(times)/(30*60):.2f} Hz')



