import sys
import json
import argparse
import numpy as np
from pathlib import Path

sys.path.append('../decode_rf_data')
from rf_data import tile_names

sys.path.append('../sat_channels')
from sat_channels import time_tree


if __name__=='__main__':

    parser = argparse.ArgumentParser(description="""
        Determine total integration time, per pointing, for all tiles
        """)

    parser.add_argument('--data_dir', metavar='\b', help='Output directoryi')
    parser.add_argument(
            '--out_dir', metavar='\b', default='./../../outputs/beam_pointings/',
            help='Output directory. Default=./../../outputs/beam_pointings/')
    
    args = parser.parse_args()
   
    data_dir    = Path(args.data_dir)
    out_dir     = Path(args.out_dir)

    # Tile names
    tiles  = tile_names()

    # lists of obs at each pointings
    with open(f'{out_dir}/obs_pointings.json') as table:
        data = json.load(table)
        start_date  = data['start_date']
        stop_date   = data['stop_date']
        point_0     = data['point_0']
        point_2     = data['point_2']
        point_4     = data['point_4']


    for tile in tiles:
        # dates: list of days
        # date_time = list of 30 min observation windows
        dates, date_time = time_tree(start_date, stop_date)
        for day in range(len(dates)):
            for window in range(len(date_time[day])):
                f_name = f'{tile}_{date_time[day][window]}.txt'
                f_path = Path(f'{data_dir}/{tile}/{dates[day]}/{f_name}')
                if f_path.is_file():
                    if date_time[day][window] in point_4:
                        print(f_path)

        break





