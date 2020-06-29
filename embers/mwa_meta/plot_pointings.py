import json
import matplotlib
import numpy as np

matplotlib.use("Agg")
import seaborn as sns
import matplotlib.pyplot as plt




if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="""
            Reads ultimate_pointing_times.json, and plots
            a histogram of total integration time at various
            pointings.
            """
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./../../outputs/beam_pointings/",
        help="Directory where json metadata lives. Default=./../../outputs/beam_pointings/",
    )
    parser.add_argument(
        "--f_name",
        metavar="\b",
        default="ultimate_pointing_times.json",
        help="File name of json to be plotted. Default=ultimate_pointing_times.json",
    )
    parser.add_argument(
        "--int_thresh",
        metavar="\b",
        default=200,
        type=int,
        help="Pointing integration time threshold. Default=200",
    )

    args = parser.parse_args()
    out_dir = args.out_dir
    f_name = args.f_name
    int_thresh = args.int_thresh

    time_point, point = point_integration(f_name, out_dir, int_thresh)
    pointing_hist(time_point, point, out_dir)
