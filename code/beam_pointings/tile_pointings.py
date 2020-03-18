import sys
import json
import argparse
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

sys.path.append('../decode_rf_data')
from rf_data import tile_names

sys.path.append('../sat_channels')
from sat_channels import time_tree


def point_int(tile, start, stop, point_0, point_2, point_4):
    
    # increment this by 0.5, for every obs which is in point_list. 
    # 0.5 = 30 min
    p_0 = 0
    p_2 = 0
    p_4 = 0

    # dates: list of days
    # date_time = list of 30 min observation windows
    dates, date_time = time_tree(start_date, stop_date)
    for day in range(len(dates)):
        for window in range(len(date_time[day])):
            f_name = f'{tile}_{date_time[day][window]}.txt'
            f_path = Path(f'{data_dir}/{tile}/{dates[day]}/{f_name}')
            if f_path.is_file():
                if date_time[day][window] in point_0:
                    p_0 += 0.5
                elif date_time[day][window] in point_2:
                    p_2 += 0.5
                elif date_time[day][window] in point_4:
                    p_4 += 0.5
                else:
                    pass
    return [p_0, p_2, p_4]

def plt_pointing_hist(sub=None, int_dict=None, tile=None):

    ax = fig.add_subplot(sub)
    
    int_list = int_dict[tile]
    x = range(len(int_list))
    leg = int_list

    pal = sns.cubehelix_palette(len(int_list), start=.4, rot=-.5, dark=.4, reverse=True)
    barplot = plt.bar(x, int_list, color=sns.color_palette(pal))
    
    def autolabel(rects):
        for idx,rect in enumerate(barplot):
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., height, 
                    leg[idx], ha='center', va='bottom', rotation=0)
    
    autolabel(barplot)
    
    plt.xticks(x, ['0','2','4'])
    plt.ylabel('Hours')
    plt.xlim(-0.7,len(int_list)-0.3)
    plt.xlabel('MWA Grid Pointing Number')
    plt.title(f'{tile}: Integration at MWA Grid Pointings')



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

    tile_integration = {}
    for tile in tiles:
        p_integration = point_int(tile, start_date, stop_date, point_0, point_2, point_4)
        tile_integration[f'{tile}'] = p_integration


    plt.style.use('seaborn')
    fig = plt.figure(figsize=(8,10))
    ax = plt_pointing_hist(sub=111, int_dict=tile_integration, tile='rf0XX')
    plt.tight_layout()
    fig.savefig('test.png')






