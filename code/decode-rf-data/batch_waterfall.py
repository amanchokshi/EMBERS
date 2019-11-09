import argparse
import rf_data as rf
from pathlib import Path
from datetime import datetime, timedelta

parser = argparse.ArgumentParser(description="""
    Plots waterfall plots for all data 
    between two date ranges. Plot for all
    tiles. Saves plots in organised directory
    structure.
    """)

parser.add_argument('--data_dir', metavar='\b', help='Dir where date is saved')
parser.add_argument('--start_date', metavar='\b', help='Date from which to start plotting waterfalls')
parser.add_argument('--stop_date', metavar='\b', help='Date until which to plot waterfalls')

args = parser.parse_args()
data_dir = args.data_dir
start_date = args.start_date
stop_date = args.stop_date

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

#data_dir = Path(data_dir)
#for i in tiles:
#    print(data_dir /i)

t = datetime.strptime('2019-11-07-17:50:00', '%Y-%m-%d-%H:%M:%S')
print(t.strftime('%H:%M'))
print((t+timedelta(hours=1)).strftime('%H:%M'))
