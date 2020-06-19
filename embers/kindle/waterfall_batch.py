"""
===============
Waterfall Batch
===============
Create waterfall plots for all rf data files recorded within a date interval.
Output files are saved to ``./embers_out/condition_data/waterfalls``

"""

import logging
import argparse
import pkg_resources
from pathlib import Path
import concurrent.futures
from itertools import repeat
from embers.condition_data.rf_data import tile_names, time_tree, batch_waterfall


parser = argparse.ArgumentParser(
    description="""
    Create waterfall plots for all rf_files within a date interval
    """
)

parser.add_argument(
    "--start_date", metavar="\b", default='', help="start date in YYYY-MM-DD format"
)
parser.add_argument(
    "--stop_date", metavar="\b", default='', help="stop date in YYYY-MM-DD format"
)
parser.add_argument(
    "--data_dir", metavar="\b", default='', help="root of dir where rf data is saved"
)
parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/condition_data",
    help="Dir where colormap sample plot is saved. Default=./embers_out/condition_data",
)

args = parser.parse_args()
start_date = args.start_date
stop_date = args.stop_date
data_dir = args.data_dir
out_dir = args.out_dir

# if no input file provided, use sample package data
if data_dir == "":
    print('----------------------------------------------------------')
    print("No data dir provided, using packaged sample data")
    data_dir = pkg_resources.resource_filename(
        "embers.kindle", "data/rf_data/"
    )

if start_date == "":
    print('No start_date provided, using 2019-10-01 for sample data')
    start_date = '2019-10-01'

if stop_date == "":
    print('No stop_date provided, using 2019-10-10 for sample data')
    print('')
    print('>>> waterfall_batch --help, for more options')
    print('----------------------------------------------------------')
    stop_date = '2019-10-10'

    # Logging config
log_dir = Path(f"{out_dir}/waterfalls")
log_dir.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    filename=f"{out_dir}/waterfalls/waterfall_batch.log",
    level=logging.INFO,
    format="%(levelname)s: %(funcName)s: %(message)s",
)


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
                    repeat(out_dir),
                )

            for result in results:
                logging.info(result)


def main():
    """Execute waterfall from terminal"""

    print(f"Processing rf data files between {start_date} and {stop_date}")
    print(f"Saving waterfall plots to: ./{log_dir}")
    waterfall_batch(start_date, stop_date, data_dir, out_dir)
