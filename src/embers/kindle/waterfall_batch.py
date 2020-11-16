"""
Waterfall Batch
---------------
"""

import argparse
import logging
from pathlib import Path

from embers.rf_tools.rf_data import waterfall_batch


def main():
    """
    Create a set of waterfall plots for all rf_files within a date interval using the :func:`~embers.rf_tools.rf_data.waterfall_batch` function.

    .. code-block:: console

        $ waterfall_batch --help

    """
    _parser = argparse.ArgumentParser(
        description="""
        Create waterfall plots for all rf_files within a date interval
        """
    )

    _parser.add_argument(
        "--start_date",
        metavar="\b",
        default="2019-10-10",
        help="start date in YYYY-MM-DD format, default=2019-10-10",
    )

    _parser.add_argument(
        "--stop_date",
        metavar="\b",
        default="2019-10-10",
        help="stop date in YYYY-MM-DD format, default=2019-10-10",
    )

    _parser.add_argument(
        "--data_dir",
        metavar="\b",
        default="./tiles_data",
        help="root of dir where rf data is saved, default=tiles_data",
    )

    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/rf_tools",
        help="Dir where colormap sample plot is saved. Default=./embers_out/rf_tools",
    )

    _args = _parser.parse_args()
    _start_date = _args.start_date
    _stop_date = _args.stop_date
    _data_dir = _args.data_dir
    _out_dir = _args.out_dir

    # Logging config
    _log_dir = Path(f"{_out_dir}/waterfalls")
    _log_dir.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=f"{_out_dir}/waterfalls/waterfall_batch.log",
        level=logging.INFO,
        format="%(levelname)s: %(funcName)s: %(message)s",
    )

    print(f"Processing rf data files between {_start_date} and {_stop_date}")
    print(f"Saving waterfall plots to: ./{_log_dir}")
    waterfall_batch(_start_date, _stop_date, _data_dir, _out_dir)
