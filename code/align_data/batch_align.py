import h5py
import argparse
import align_data as ali
import concurrent.futures

from pathlib import Path
from itertools import product
from datetime import datetime, timedelta

#TODO needed to add code dir to PYTHONPATH for this import to work. 
#TODO Is this the best way?
from decode_rf_data import rf_data as rf


parser = argparse.ArgumentParser(description="""
    Time alignes and smoothes all data 
    between two date ranges. Saves data
    pairs in organised directory 
    structure as hdf5 files.
    """)

parser.add_argument('--data_dir', metavar='\b', help='Dir where date is saved')
parser.add_argument('--start_date', metavar='\b', help='Date from which to start aligning data. Ex: 2019-10-10')
parser.add_argument('--stop_date', metavar='\b', help='Date until which to align data. Ex: 2019-10-11')
parser.add_argument('--out_dir', metavar='\b', default='./../../outputs/align_data/',help='Output directory. Default=./../../outputs/align_data/')

args = parser.parse_args()
data_dir = args.data_dir
start_date = args.start_date
stop_date = args.stop_date
out_dir = args.out_dir

# Import list of tile names from rf_data.py
tiles = rf.tile_names()

refs = tiles[:4]
AUTS = tiles[4:]

# Split the list into XX and YY lists
refs_XX = refs[::2]
refs_YY = refs[1::2]

AUTS_XX = AUTS[::2]
AUTS_YY = AUTS[1::2]


# Create a list of pairs of tiles to be aligned. 
# Essentially all possible combinations of Refs, AUTS
align_pairs = []
for pair in list(product(refs_XX, AUTS_XX)):
    align_pairs.append(pair)

for pair in list(product(refs_YY, AUTS_YY)):
    align_pairs.append(pair)


# Time stuff to help traverse data dir tree.
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


def save_align_data():
    ref_name = f'{ref}_{time_stamp}'
    aut_name = f'{aut}_{time_stamp}'

#for R


#with h5py.File('test.hdf5', 'w') as f:
#    f.create_dataset('ref_p_aligned', data=ref_p_aligned)
#    f.create_dataset('tile_p_aligned', data=tile_p_aligned)
#    f.create_dataset('time_array', data=time_array)

