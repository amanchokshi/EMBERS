import json
import numpy as np
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt

import sys

sys.path.append("../sat_channels")
from window_chan_map import time_tree


def point_int(
    tile, start_date, stop_date, point_0, point_2, point_4, point_41, data_dir
):
    """Calculate integration at each pointing for every tile"""

    # increment this by 0.5, for every obs which is in point_list.
    # 0.5 = 30 min
    p_0 = 0
    p_2 = 0
    p_4 = 0
    p_41 = 0

    # dates: list of days
    # date_time = list of 30 min observation windows
    dates, date_time = time_tree(start_date, stop_date)
    for day in range(len(dates)):
        for window in range(len(date_time[day])):
            f_name = f"{tile}_{date_time[day][window]}.txt"
            f_path = Path(f"{data_dir}/{tile}/{dates[day]}/{f_name}")
            if f_path.is_file():
                if date_time[day][window] in point_0:
                    p_0 += 0.5
                elif date_time[day][window] in point_2:
                    p_2 += 0.5
                elif date_time[day][window] in point_4:
                    p_4 += 0.5
                elif date_time[day][window] in point_41:
                    p_41 += 0.5
                else:
                    pass
    return [p_0, p_2, p_4, p_41]


def plt_pointing_hist(a, b, c, fig, int_dict=None, tile=None, x_ticks=None):

    ax = fig.add_subplot(a, b, c)

    y_max = max(i for v in int_dict.values() for i in v)
    int_list = int_dict[tile]
    x = range(len(int_list))
    leg = int_list

    pal = sns.cubehelix_palette(
        len(int_list), start=0.7, rot=-0.7, dark=0.4, reverse=True
    )
    barplot = plt.bar(x, int_list, color=sns.color_palette(pal), edgecolor="black")

    def autolabel(rects):
        for idx, rect in enumerate(barplot):
            height = rect.get_height()
            ax.text(
                rect.get_x() + rect.get_width() / 2.0,
                height,
                leg[idx],
                ha="center",
                va="bottom",
                rotation=0,
            )

    autolabel(barplot)

    # these are matplotlib.patch.Patch properties
    props = dict(boxstyle="round", facecolor="wheat", alpha=0.5)

    # place a text box in upper left in axes coords
    ax.text(
        0.65,
        0.95,
        f"{tile}",
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=props,
    )

    plt.xlim(-0.5, len(int_list) - 0.5)
    ax.set_yticklabels([])

    ax.set_xticklabels([])
    ax.set_ylim(0, 1.1 * y_max)

    return ax


if __name__ == "__main__":

    import sys
    import argparse

    sys.path.append("../decode_rf_data")
    from rf_data import tile_names

    parser = argparse.ArgumentParser(
        description="""
        Determine total integration time, per pointing, for all tiles
        """
    )

    parser.add_argument("--data_dir", metavar="\b", help="Data directory")
    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./../../outputs/beam_pointings/",
        help="Output directory. Default=./../../outputs/beam_pointings/",
    )

    args = parser.parse_args()

    data_dir = Path(args.data_dir)
    out_dir = Path(args.out_dir)

    # Tile names
    tiles = tile_names()

    # lists of obs at each pointings
    with open(f"{out_dir}/obs_pointings.json") as table:
        data = json.load(table)
        start_date = data["start_date"]
        stop_date = data["stop_date"]
        point_0 = data["point_0"]
        point_2 = data["point_2"]
        point_4 = data["point_4"]
        point_41 = data["point_41"]

    tile_integration = {}
    for tile in tiles:
        p_integration = point_int(
            tile, start_date, stop_date, point_0, point_2, point_4, point_41, data_dir
        )
        tile_integration[f"{tile}"] = p_integration

    plt.style.use("seaborn")
    fig = plt.figure(figsize=(14, 10))

    ax01 = plt_pointing_hist(4, 8, 1, fig, int_dict=tile_integration, tile=tiles[0])
    ax02 = plt_pointing_hist(4, 8, 2, fig, int_dict=tile_integration, tile=tiles[4])
    ax03 = plt_pointing_hist(4, 8, 3, fig, int_dict=tile_integration, tile=tiles[5])
    ax04 = plt_pointing_hist(4, 8, 4, fig, int_dict=tile_integration, tile=tiles[6])
    ax05 = plt_pointing_hist(4, 8, 5, fig, int_dict=tile_integration, tile=tiles[7])
    ax06 = plt_pointing_hist(4, 8, 6, fig, int_dict=tile_integration, tile=tiles[8])
    ax07 = plt_pointing_hist(4, 8, 7, fig, int_dict=tile_integration, tile=tiles[9])
    ax08 = plt_pointing_hist(4, 8, 8, fig, int_dict=tile_integration, tile=tiles[10])
    ax09 = plt_pointing_hist(4, 8, 9, fig, int_dict=tile_integration, tile=tiles[1])
    ax10 = plt_pointing_hist(4, 8, 10, fig, int_dict=tile_integration, tile=tiles[11])
    ax11 = plt_pointing_hist(4, 8, 11, fig, int_dict=tile_integration, tile=tiles[12])
    ax12 = plt_pointing_hist(4, 8, 12, fig, int_dict=tile_integration, tile=tiles[13])
    ax12 = plt_pointing_hist(4, 8, 13, fig, int_dict=tile_integration, tile=tiles[14])
    ax13 = plt_pointing_hist(4, 8, 14, fig, int_dict=tile_integration, tile=tiles[15])
    ax14 = plt_pointing_hist(4, 8, 15, fig, int_dict=tile_integration, tile=tiles[16])
    ax15 = plt_pointing_hist(4, 8, 16, fig, int_dict=tile_integration, tile=tiles[17])
    ax16 = plt_pointing_hist(4, 8, 17, fig, int_dict=tile_integration, tile=tiles[2])
    ax17 = plt_pointing_hist(4, 8, 18, fig, int_dict=tile_integration, tile=tiles[18])
    ax18 = plt_pointing_hist(4, 8, 19, fig, int_dict=tile_integration, tile=tiles[19])
    ax19 = plt_pointing_hist(4, 8, 20, fig, int_dict=tile_integration, tile=tiles[20])
    ax20 = plt_pointing_hist(4, 8, 21, fig, int_dict=tile_integration, tile=tiles[21])
    ax21 = plt_pointing_hist(4, 8, 22, fig, int_dict=tile_integration, tile=tiles[22])
    ax22 = plt_pointing_hist(4, 8, 23, fig, int_dict=tile_integration, tile=tiles[23])
    ax23 = plt_pointing_hist(4, 8, 24, fig, int_dict=tile_integration, tile=tiles[24])
    ax24 = plt_pointing_hist(4, 8, 25, fig, int_dict=tile_integration, tile=tiles[3])
    ax25 = plt_pointing_hist(4, 8, 26, fig, int_dict=tile_integration, tile=tiles[25])
    ax26 = plt_pointing_hist(4, 8, 27, fig, int_dict=tile_integration, tile=tiles[26])
    ax27 = plt_pointing_hist(4, 8, 28, fig, int_dict=tile_integration, tile=tiles[27])
    ax28 = plt_pointing_hist(4, 8, 29, fig, int_dict=tile_integration, tile=tiles[28])
    ax29 = plt_pointing_hist(4, 8, 30, fig, int_dict=tile_integration, tile=tiles[29])
    ax30 = plt_pointing_hist(4, 8, 31, fig, int_dict=tile_integration, tile=tiles[30])
    ax31 = plt_pointing_hist(4, 8, 32, fig, int_dict=tile_integration, tile=tiles[31])

    #    plt.title('Pointing Integration per Tile')
    #    plt.ylabel('Hours')
    #    plt.xlabel('Pointing [0, 2, 4]')
    #    fig.text(0.5, 0.04, 'MWA Receiver Pointing [0, 2, 4]', ha='center')
    #    fig.text(0.04, 0.5, 'Hours', va='center', rotation='vertical')
    plt.tight_layout()
    fig.savefig(f"{out_dir}/tiles_pointing_integration.png")
