import os
import sys
import rf_data as rf
from pathlib import Path
import concurrent.futures
from itertools import repeat
from datetime import datetime, timedelta


def plt_single_waterfall(rf_dir, rf_filename, out_dir):
    """Plot a single waterfall"""

    power, times = rf.read_data(f'{rf_dir}/{rf_filename}.txt')
    plt = rf.plot_waterfall(power, times, rf_filename)
    
    # Make output dir if it doesn't exist
    os.makedirs(os.path.dirname(out_dir), exist_ok=True)
    
    plt.savefig(f'{out_dir}/{rf_filename}.png')
    plt.close()


def waterfall_plot(tile, out_dir, data_dir, date, time_stamp):
    rf_name = f'{tile}_{time_stamp}'
    rf_path = Path(f'{data_dir}/{tile}/{date}')
    
    try:
        power, times = rf.read_data(f'{rf_path}/{rf_name}.txt')

        plt = rf.plot_waterfall(power, times, rf_name)

        save_dir = Path(f'{out_dir}/waterfalls/{date}/{time_stamp}')
        save_dir.mkdir(parents=True, exist_ok=True)
        
        
        plt.savefig(f'{save_dir}/{rf_name}.png')
        plt.close()
        return f'Waterfall {rf_name}.png saved'
    except Exception:
        return f'File {rf_name}.txt missing'

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="""
        Plots waterfall plots for all data 
        between two date ranges. Plot for all
        tiles. Saves plots in organised directory
        structure.
        """)
    
    parser.add_argument('--data_dir', metavar='\b', help='Dir where date is saved')
    parser.add_argument('--start_date', metavar='\b', help='Date from which to start plotting waterfalls. Ex: 2019-10-10')
    parser.add_argument('--stop_date', metavar='\b', help='Date until which to plot waterfalls. Ex: 2019-10-11')
    parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/decode_rf_data',help='Output directory. Default=./../../outputs/decode_rf_data/')
    
    args = parser.parse_args()
    data_dir = args.data_dir
    start_date = args.start_date
    stop_date = args.stop_date
    out_dir = args.out_dir
    
    # Save logs 
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f'{out_dir}/logs_{start_date}_{stop_date}.txt', 'a')
    
    # Import list of tile names from rf_data.py
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
            
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(
                        waterfall_plot,
                        repeat(tile),
                        repeat(out_dir), 
                        repeat(data_dir), 
                        repeat(dates[d]), 
                        date_time[d])
               
            for result in results:
                print(result)



