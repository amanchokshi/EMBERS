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


_parser = argparse.ArgumentParser(
    description="""
    Saved a single waterfall plot
    """
)

_parser.add_argument(
    "--rf_file", metavar="\b", default="", help="Path to raw rf data file"
)
_parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="embers_out/condition_data",
    help="Dir where colormap sample plot is saved. Default=./embers_out/condition_data",
)

_args = _parser.parse_args()
_rf_file = _args.rf_file
_out_dir = _args.out_dir

# if no input file provided, use sample package data
if _rf_file == "":
    print("--------------------------------------------------")
    print("No input data provided, using packaged sample data")
    print(">>> waterfall --help, for more options")
    print("--------------------------------------------------")
    _rf_file = pkg_resources.resource_filename(
        "embers.kindle", "data/rf_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt"
    )


def main():
    """Execute waterfall from terminal."""

    _rf_name = Path(_rf_file).stem
    print(f"Waterfall plot saved to ./{_out_dir}/{_rf_name}.png")
    single_waterfall(_rf_file, _out_dir)
