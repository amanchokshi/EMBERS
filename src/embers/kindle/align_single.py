"""
Align Single
============

Plot a single frequency channel of raw rf data and the result of
interpolation and Savitzky–Golay smoothing. Save plot to :samp:`./embers_out/rf_tools`

"""

import argparse

from embers.rf_tools.align_data import plot_savgol_interp


def main():
    """
    Temporally align reference and tile data using the :func:`~embers.rf_tools.align_data.plot_savgol_interp` function.

    .. code-block:: console

        $ align_single --help

    """
    _parser = argparse.ArgumentParser(
        description="""
        Create sample savgol interp plot of pair of rf data
        """
    )

    _parser.add_argument(
        "--ref_file",
        metavar="\b",
        default="tiles_data/rf0XX/2019-10-10/rf0XX_2019-10-10-02:30.txt",
        help="Reference rf data file. Default=tiles_data/rf0XX/2019-10-10/rf0XX_2019-10-10-02:30.txt",
    )
    _parser.add_argument(
        "--tile_file",
        metavar="\b",
        default="tiles_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt",
        help="Tile rf data file. Default=tiles_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt",
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
        "--channel",
        metavar="\b",
        default=59,
        type=int,
        help="Frequency channel to plot",
    )
    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/rf_tools",
        help="Dir where savgol-interp sample plot is saved. Default=./embers_out/rf_tools",
    )

    _args = _parser.parse_args()
    _ref_file = _args.ref_file
    _tile_file = _args.tile_file
    _savgol_window_1 = _args.savgol_window_1
    _savgol_window_2 = _args.savgol_window_2
    _polyorder = _args.polyorder
    _interp_type = _args.interp_type
    _interp_freq = _args.interp_freq
    _channel = _args.channel
    _out_dir = _args.out_dir

    print(f"Saving sample savgol_interp plot to: {_out_dir}/savgol_interp_sample.png")
    plot_savgol_interp(
        ref=_ref_file,
        tile=_tile_file,
        savgol_window_1=_savgol_window_1,
        savgol_window_2=_savgol_window_2,
        polyorder=_polyorder,
        interp_type=_interp_type,
        interp_freq=_interp_freq,
        channel=_channel,
        out_dir=_out_dir,
    )
