"""
=========
Waterfall
=========
Saves a single waterfall plot to ``./embers_out/condition_data``

"""

import argparse
from pathlib import Path
from embers.condition_data.rf_data import single_waterfall

parser = argparse.ArgumentParser(
    description="""
    Saved a single waterfall plot
    """
)

parser.add_argument(
    "--rf_file", metavar="\b", required=True, help="Path to raw rf data file"
)
parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/condition_data",
    help="Dir where colormap sample plot is saved. Default=./embers_out/condition_data",
)

args = parser.parse_args()
rf_file = args.rf_file
out_dir = args.out_dir


def main():
    """Execute waterfall from terminal"""

    rf_name = Path(rf_file).stem
    print(f"Waterfall plot saved to ./{out_dir}/{rf_name}.png")
    single_waterfall(rf_file, out_dir)
