"""
===============
Waterfall Batch
===============
Create waterfall plots for all rf data files recorded within a date interval.
Output files are saved to ``./embers_out/condition_data/waterfalls``

"""

import argparse
from pathlib import Path
import concurrent.futures
from itertools import repeat
from embers.condition_data.rf_data import tile_names, time_tree, batch_waterfall

parser = argparse.ArgumentParser(
    description="""
    Create waterfall plots for all rf_files within a date interval
    """
)

parser.add_argument("--start_date", metavar="\b", help="start date in YYYY-MM-DD format")
parser.add_argument("--stop_date", metavar="\b", help="stop date in YYYY-MM-DD format")
parser.add_argument(
    "--data_dir",
    metavar="\b",
    help="root of dir where rf data is saved"
)
parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/condition_data",
    help="Dir where colormap sample plot is saved. Default=./embers_out/condition_data"
)

args = parser.parse_args()
start_date = args.start_date
stop_date = args.stop_date
data_dir = args.data_dir
out_dir = args.out_dir



def waterfall_batch(start_date, stop_date, data_dir, out_dir):

    dates, time_stamps = time_tree(start_date, stop_date)
    
    for tile in tile_names():
        for day in range(len(dates)):
    
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(
                    batch_waterfall,
                    repeat(tile),
                    time_stamps[day],
                    repeat(data_dir),
                    repeat(out_dir)
                )
    
            for result in results:
                if result != 'cannot unpack non-iterable NoneType object':
                    print(result)


waterfall_batch(start_date, stop_date, data_dir, out_dir)
#def main():
#    """Execute waterfall from terminal"""
#
#    rf_name = Path(rf_file).stem
#    print(f"Waterfall plot saved to ./{out_dir}/{rf_name}.png")
#    single_waterfall(rf_file, out_dir)
