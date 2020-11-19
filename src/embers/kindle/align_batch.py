"""
Align Batch
===========

Interpolate and Savitzky–Golay smooth a pair of raw rf data files
before temorally aligning the. Save aligned data as npz files to
``./embers_out/rf_tools/align_data``

"""

import argparse

from embers.rf_tools.align_data import align_batch


def main():
    """
    Temporally align all RF files within a date interval using the :func:`~embers.rf_tools.align_data.align_batch` function.

    .. code-block:: console

        $ align_batch --help
    """

    _parser = argparse.ArgumentParser(
        description="""
        Batch align pairs of rf data files between a date interval
        """
    )

    _parser.add_argument(
        "--start_date",
        metavar="\b",
        default="2019-10-10",
        help="start date in YYYY-MM-DD format",
    )

    _parser.add_argument(
        "--stop_date",
        metavar="\b",
        default="2019-10-10",
        help="stop date in YYYY-MM-DD format",
    )

    _parser.add_argument(
        "--savgol_window_1",
        metavar="\b",
        type=int,
        default=11,
        help="First savgol window. Default=11",
    )

    _parser.add_argument(
        "--savgol_window_2",
        metavar="\b",
        type=int,
        default=15,
        help="Second savgol window. Default=15",
    )

    _parser.add_argument(
        "--polyorder",
        metavar="\b",
        type=int,
        default=2,
        help="Polynomial order. Default=2",
    )

    _parser.add_argument(
        "--interp_type",
        metavar="\b",
        default="cubic",
        help="Interpolation type. Default=cubic",
    )

    _parser.add_argument(
        "--interp_freq",
        metavar="\b",
        type=int,
        default=1,
        help="Interpolation frequency. Default=1",
    )

    _parser.add_argument(
        "--data_dir",
        metavar="\b",
        default="./tiles_data",
        help="Root of dir where rf data is saved. Default=./tiles_data",
    )

    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/rf_tools/align_data",
        help="Dir where savgol-interp sample plot is saved. Default=./embers_out/rf_tools/align_data",
    )

    _parser.add_argument(
        "--max_cores",
        metavar="\b",
        type=int,
        help="Maximum number of cores to be used by this script. By default all core available cores are used",
    )

    _args = _parser.parse_args()
    _start_date = _args.start_date
    _stop_date = _args.stop_date
    _savgol_window_1 = _args.savgol_window_1
    _savgol_window_2 = _args.savgol_window_2
    _polyorder = _args.polyorder
    _interp_type = _args.interp_type
    _interp_freq = _args.interp_freq
    _data_dir = _args.data_dir
    _out_dir = _args.out_dir
    _max_cores = _args.max_cores

    print(f"Aligned files saved to: {_out_dir}")
    align_batch(
        start_date=_start_date,
        stop_date=_stop_date,
        savgol_window_1=_savgol_window_1,
        savgol_window_2=_savgol_window_2,
        polyorder=_polyorder,
        interp_type=_interp_type,
        interp_freq=_interp_freq,
        data_dir=_data_dir,
        out_dir=_out_dir,
        max_cores=_max_cores,
    )
