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
            '--out_dir', metavar='\b', default='./../../outputs/tile_maps/',
            help='Output directory. Default=./../../outputs/tile_maps/')
    
    args = parser.parse_args()
   
    data_dir    = Path(args.data_dir)
    out_dir     = Path(args.out_dir)

    # Tile names
    tiles  = tile_names()

    # lists of obs at each pointings
    with open(f'{out_dir}/{f_name}') as table:
        data = json.load(table)
        point_0 = data['point_0']
        point_2 = data['point_2']
        point_4 = data['point_4']

    # dates: list of days
    # date_time = list of 30 min observation windows
    dates, date_time = time_tree(start_date, stop_date)
    for day in range(len(dates)):
        for window in range(len(date_time[day])):
            for x in tilesx:
                if Path(f'{align_dir}/{dates[day]}/{date_time[day][window]}/rf0XX_{x}_{date_time[day][window]}_aligned.npz').is_file():
                    print(f'rf0XX_{x}_{date_time[day][window]}_aligned.npz exists')
                else:
                    print(f'rf0XX_{x}_{date_time[day][window]}_aligned.npz missing')
            #print(date_time[day][window])





