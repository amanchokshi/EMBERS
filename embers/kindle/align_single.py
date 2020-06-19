"""
============
Align Single
============

Plot a single frequency channel of raw rf data and the result of
interpolation and Savitzky–Golay smoothing. Save plot to ``./embers_out/rf_tools``

"""

import argparse
import pkg_resources
from pathlib import Path
import matplotlib.pyplot as plt
from embers.rf_tools.align_data import savgol_interp

_parser = argparse.ArgumentParser(
    description="""
    Create sample savgol interp plot of pair of rf data
    """
)

_parser.add_argument(
    "--ref_file", metavar="\b", default="", help="Reference rf data file"
)

_parser.add_argument("--tile_file", metavar="\b", default="", help="Tile rf data file")

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
    "--polyorder", metavar="\b", type=int, default=2, help="Polynomial order. Default=2"
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
    "--channel", metavar="\b", default="", help="Frequency channel to plot"
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


# if no input files provided, use sample package data
if _ref_file == "":
    print("----------------------------------------------------------")
    print("No ref_file provided, using packaged sample data")
    _ref_file = pkg_resources.resource_filename(
        "embers.kindle", "data/rf_data/rf0XX/2019-10-10/rf0XX_2019-10-10-02:30.txt"
    )

if _tile_file == "":
    print("No tile_file provided, using packaged sample data")
    _tile_file = pkg_resources.resource_filename(
        "embers.kindle", "data/rf_data/S06XX/2019-10-10/S06XX_2019-10-10-02:30.txt"
    )

if _channel == "":
    print("No frequency channel provided, using 59 for sample data")
    print("")
    print(">>> savgol_interp_sample --help, for more options")
    print("----------------------------------------------------------")
    _channel = 59


def plot_savgol_interp(
    ref=None,
    tile=None,
    savgol_window_1=None,
    savgol_window_2=None,
    polyorder=None,
    interp_type=None,
    interp_freq=None,
    channel=None,
    out_dir=None,
):

    (
        ref_ali,
        tile_ali,
        time_array,
        ref_power,
        tile_power,
        ref_time,
        tile_time,
    ) = savgol_interp(
        ref=ref,
        tile=tile,
        savgol_window_1=savgol_window_1,
        savgol_window_2=savgol_window_2,
        polyorder=polyorder,
        interp_type=interp_type,
        interp_freq=interp_freq,
    )

    # Sample align plot
    plt.style.use("seaborn")
    plt.rcParams["figure.figsize"] = (9, 6)

    # convert times to minuts from first datapoint
    time_array = (time_array - time_array[0]) / 60
    ref_time = (ref_time - ref_time[0]) / 60
    tile_time = (tile_time - tile_time[0]) / 60

    plt.plot(
        time_array,
        tile_ali[::, channel],
        color="#e23a4e",
        alpha=0.9,
        label="tile savgol",
    )
    plt.scatter(
        tile_time,
        tile_power[::, channel],
        color="#f78b51",
        marker=".",
        alpha=0.6,
        label="tile raw",
    )

    plt.plot(
        time_array,
        ref_ali[::, channel],
        color="#252b40",
        alpha=0.9,
        label="ref savgol",
    )
    plt.scatter(
        ref_time,
        ref_power[::, channel],
        color="#6a82bb",
        marker=".",
        alpha=0.6,
        label="ref raw",
    )

    leg = plt.legend(loc="upper left", frameon=True)
    leg.get_frame().set_facecolor("white")
    for l in leg.legendHandles:
        l.set_alpha(1)

    plt.ylim(-110, -20)
    plt.ylabel("Raw Power [dBm]")
    plt.xlabel("Time [min]")
    plt.tight_layout()
    Path(f"{out_dir}").mkdir(parents=True, exist_ok=True)
    plt.savefig(f"{out_dir}/savgol_interp_sample.png")


def main():
    """Execute align_single from terminal."""

    print(f"Saving sample savgol_interp plot to: {_out_dir}")
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
