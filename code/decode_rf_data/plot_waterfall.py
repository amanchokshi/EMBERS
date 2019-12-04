import os
import argparse
import rf_data as rf


parser = argparse.ArgumentParser(description="""
        Will plot a single waterfall plot
        """)

parser.add_argument('--rf_dir', metavar='\b', default='./../../data/', help='Path to rf data directory. Default=./../../data/')
parser.add_argument('--rf_name', metavar='\b', default='rf0XX_2019-10-10-02:30', help='Name of rf data file. Default=rf0XX_2019-10-10-02:30')
parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/decode_rf_data/', help='Output dir. Default=./../../outputs/decode_rf_dir/')

args = parser.parse_args()
rf_dir = args.rf_dir
rf_name = args.rf_name
out_dir = args.out_dir


power, times = rf.read_data(f'{rf_dir}/{rf_name}.txt')
plt = rf.plot_waterfall(power, times, rf_name)

# Make output dir if it doesn't exist
os.makedirs(os.path.dirname(out_dir), exist_ok=True)

plt.savefig('{out_dir}/{rf_name}.png')

