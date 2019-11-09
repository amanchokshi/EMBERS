import argparse
import rf_data as rf


parser = argparse.ArgumentParser(description="""
        Will plot a single waterfall plot
        """)

parser.add_argument('--rf_dir', metavar='\b', default='./../../sample-data/', help='Path to rf data directory. Default=./../../sample-data/')
parser.add_argument('--rf_name', metavar='\b', default='rf0XX_2019-10-10-15:00', help='Name of rf data file. Default=rf0XX_2019-10-10-15:00')

args = parser.parse_args()
rf_dir = args.rf_dir
rf_name = args.rf_name


power, times = rf.read_data('{}{}.txt'.format(rf_dir,rf_name))
plt = rf.plot_waterfall(power, times, rf_name)
plt.show()


