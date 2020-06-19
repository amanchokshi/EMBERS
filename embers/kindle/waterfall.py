"""
=========
Waterfall
=========
Saves a single waterfall plot to ``./embers_out/condition_data``

"""

import argparse
import pkg_resources
from pathlib import Path
from embers.condition_data.rf_data import single_waterfall


parser = argparse.ArgumentParser(
    description="""
    Saved a single waterfall plot
    """
)

parser.add_argument(
    "--rf_file", metavar="\b", default="", help="Path to raw rf data file"
)
parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="embers_out/condition_data",
    help="Dir where colormap sample plot is saved. Default=./embers_out/condition_data",
)

args = parser.parse_args()
rf_file = args.rf_file
out_dir = args.out_dir

# if no input file provided, use sample package data
if rf_file == "":
    print("No input data provided, using packaged sample data")
    rf_file = pkg_resources.resource_filename(
        "embers.kindle", "data/rf_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt"
    )


def main():
    """Execute waterfall from terminal"""

    rf_name = Path(rf_file).stem
    print(f"Waterfall plot saved to ./{out_dir}/{rf_name}.png")
    single_waterfall(rf_file, out_dir)
