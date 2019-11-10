import argparse
import rf_data as rf
from pathlib import Path
import concurrent.futures
from datetime import datetime, timedelta

parser = argparse.ArgumentParser(description="""
    Plots waterfall plots for all data 
    between two date ranges. Plot for all
    tiles. Saves plots in organised directory
    structure.
    """)

parser.add_argument('--data_dir', metavar='\b', help='Dir where date is saved')
parser.add_argument('--start_date', metavar='\b', help='Date from which to start plotting waterfalls. Ex: 2019-10-10')
parser.add_argument('--stop_date', metavar='\b', help='Date until which to plot waterfalls. Ex: 2019-10-11')
parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/decode-rf-data',help='Output directory. Default=./../../outputs/decode-rf-data/')

args = parser.parse_args()
data_dir = args.data_dir
start_date = args.start_date
stop_date = args.stop_date
out_dir = args.out_dir

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

#for tile in tiles:
#    for d in range(len(dates)):
#        rf_path = data_dir/tile/dates[d]
#        for time_stamp in date_time[d]:
#            try:
#                rf_name = f'{tile}_{time_stamp}'
#                power, times = rf.read_data(f'{rf_path}/{rf_name}.txt')
#                plt = rf.plot_waterfall(power, times, rf_name)
#
#                save_dir = out_dir/'waterfalls'/dates[d]/time_stamp
#                save_dir.mkdir(parents=True, exist_ok=True)
#                
#                
#                plt.savefig(f'{save_dir}/{rf_name}.png')
#                plt.close()
#            except Exception:
#                print(f'File {rf_name}.txt missing')
           


def waterfall_plot(time_stamp):
    try:
        rf_name = f'{tile}_{time_stamp}'
        power, times = rf.read_data(f'{rf_path}/{rf_name}.txt')
        plt = rf.plot_waterfall(power, times, rf_name)

        save_dir = out_dir/'waterfalls'/dates[d]/time_stamp
        save_dir.mkdir(parents=True, exist_ok=True)
        
        
        plt.savefig(f'{save_dir}/{rf_name}.png')
        plt.close()
    except Exception:
        return f'File {rf_name}.txt missing'

for tile in tiles:
    for d in range(len(dates)):
        rf_path = data_dir/tile/dates[d]
        
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(waterfall_plot, date_time[d])
           
        for result in results:
            print(result)

#        for time_stamp in date_time[d]:
#            waterfall_plot(time_stamp)



