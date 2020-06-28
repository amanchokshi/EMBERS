"""
Satellite Channels
------------------

Create satellite channel maps which identify the transmission 
channel [frequency] of every satellite in a set of 30 minute observations.

"""

import sys
import argparse
from pathlib import Path
from embers.sat_utils.sat_channels import batch_window_map


_parser = argparse.ArgumentParser(
    description="""
    Create a batch of window channel maps for rf observations between a date interval
    """
)
_parser.add_argument(
    "--start_date", metavar="\b", default="", help="start date in YYYY-MM-DD format"
)
_parser.add_argument(
    "--stop_date", metavar="\b", default="", help="stop date in YYYY-MM-DD format"
)
_parser.add_argument(
    "--ali_dir",
    metavar="\b",
    default="./embers_out/rf_tools/align_data",
    help="root of dir where aligned npz files are saved. Default:./embers_out/rf_tools/align_data",
)
_parser.add_argument(
    "--chrono_dir",
    metavar="\b",
    default="./embers_out/sat_utils/ephem_chrono",
    help="root of dir where chrono ephem json files are saved. Default:./embers_out/sat_utils/ephem_chrono",
)
_parser.add_argument(
    "--sat_thresh",
    metavar="\b",
    type=int,
    default=1,
    help="Satellite threshold. Integer multiple of standard deviations of rf power array, used to exclude occupied channels when determining noise floor. Default:1",
)
_parser.add_argument(
    "--noi_thresh",
    metavar="\b",
    type=int,
    default=3,
    help="Noise threshold. Integer multiple of MAD of noisy data, used to determining noise floor. Default:3",
)
_parser.add_argument(
    "--pow_thresh",
    metavar="\b",
    type=int,
    default=15,
    help="Power threshold, used to classify channels with satellites. dB above noise floor. Default:15",
)
_parser.add_argument(
    "--occ_thresh",
    metavar="\b",
    type=float,
    default=0.8,
    help="Minimum window occupation to classify channels with satellites. Default:0.80",
)
_parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/sat_utils/sat_channels",
    help="root of dir where window channel maps will be saved. Default:./embers_out/sat_utils/sat_channels",
)

_parser.add_argument(
    "--plots",
    metavar="\b",
    default="False",
    help="If True, create a bunch of useful diagnostic plots. Default=False",
)

_args = _parser.parse_args()
_start_date = _args.start_date
_stop_date = _args.stop_date
_ali_dir = _args.ali_dir
_chrono_dir = _args.chrono_dir
_sat_thresh = _args.sat_thresh
_noi_thresh = _args.noi_thresh
_pow_thresh = _args.pow_thresh
_occ_thresh = _args.occ_thresh
_out_dir = _args.out_dir
_plots = _args.plots

if _plots == "True":
    _plots = True

Path(_out_dir).mkdir(parents=True, exist_ok=True)


# if no input file provided, use sample package data
if not Path(_ali_dir).is_dir():
    print("-------------------------------------------")
    print("Check input data directories")
    print(">>> sat_channels --help, for more options")
    print("-------------------------------------------")


def main():
    """Execute sat_channels from terminal."""

    print(f"Window channel maps will be saved to: {_out_dir}")
    sys.stdout = open(f"{_out_dir}/sat_channels.log", "a")

    batch_window_map(
        _start_date,
        _stop_date,
        _ali_dir,
        _chrono_dir,
        _sat_thresh,
        _noi_thresh,
        _pow_thresh,
        _occ_thresh,
        _out_dir,
        plots=_plots,
    )
