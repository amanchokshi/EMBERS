import numpy as np
from align_data import savgol_interp
import concurrent.futures

from pathlib import Path
from itertools import repeat
from itertools import product
from datetime import datetime, timedelta

import sys

sys.path.append("../decode_rf_data")
import rf_data as rf


def save_aligned(
    ref,
    aut,
    ref_path,
    aut_path,
    savgol_window_1,
    savgol_window_2,
    polyorder,
    interp_type,
    interp_freq,
    out_dir,
    date,
    time_stamp,
):

    try:
        ref_file = f"{ref_path}/{ref}_{time_stamp}.txt"
        aut_file = f"{aut_path}/{aut}_{time_stamp}.txt"

        _, _, _, _, ref_p_aligned, tile_p_aligned, time_array = savgol_interp(
            ref_file,
            aut_file,
            savgol_window_1=savgol_window_1,
            savgol_window_2=savgol_window_2,
            polyorder=polyorder,
            interp_type=interp_type,
            interp_freq=interp_freq,
        )

        # creates output directory if it doesn't exist
        save_dir = Path(f"{out_dir}/{date}/{time_stamp}")
        save_dir.mkdir(parents=True, exist_ok=True)

        # Convert the power array to float32
        # Convert list of times to float64 (double)
        # Save as compressed npz file. Seems to drastically reduce size
        np.savez_compressed(
            f"{save_dir}/{ref}_{aut}_{time_stamp}_aligned.npz",
            ref_p_aligned=np.single(ref_p_aligned),
            tile_p_aligned=np.single(tile_p_aligned),
            time_array=np.double(time_array),
        )

        return f"Saving {ref}_{aut}_{time_stamp}_aligned.npz"

    except Exception:
        return f"Cound not save {ref}_{aut}_{time_stamp}_aligned.npz. Missing file"
        # This exception should only be raised if both files don't exist
        # In that case, the exception from rf_data.read_data() will be
        # displayed here.


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="""
        Time alignes and smoothes all data 
        between two date ranges. Saves data
        pairs in organised directory 
        structure as .npz files.
        """
    )

    parser.add_argument("--data_dir", metavar="\b", help="Dir where date is saved")
    parser.add_argument(
        "--start_date",
        metavar="\b",
        help="Date from which to start aligning data. Ex: 2019-10-10",
    )
    parser.add_argument(
        "--stop_date",
        metavar="\b",
        help="Date until which to align data. Ex: 2019-10-11",
    )
    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./../../outputs/align_data/",
        help="Output directory. Default=./../../outputs/align_data/",
    )
    parser.add_argument(
        "--savgol_window_1",
        metavar="\b",
        type=int,
        default=11,
        help="Length of savgol window. Must be odd. Default=11",
    )
    parser.add_argument(
        "--savgol_window_2",
        metavar="\b",
        type=int,
        default=15,
        help="Length of savgol window. Must be odd. Default=15",
    )
    parser.add_argument(
        "--polyorder",
        metavar="\b",
        type=int,
        default=2,
        help="Order of polynomial to fit to savgol window. Default=2",
    )
    parser.add_argument(
        "--interp_type",
        metavar="\b",
        default="cubic",
        help="Type of interpolation. Ex: cubic, linear, etc. Default=cubic",
    )
    parser.add_argument(
        "--interp_freq",
        metavar="\b",
        type=int,
        default=1,
        help="Frequency at which to resample smoothed data, in Hertz. Default=1",
    )

    args = parser.parse_args()

    data_dir = args.data_dir
    start_date = args.start_date
    stop_date = args.stop_date
    out_dir = args.out_dir
    savgol_window_1 = args.savgol_window_1
    savgol_window_2 = args.savgol_window_2
    polyorder = args.polyorder
    interp_type = args.interp_type
    interp_freq = args.interp_freq

    # Save logs
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f"{out_dir}/logs_{start_date}_{stop_date}.txt", "a")

    # Import list of tile names from rf_data.py
    tiles = rf.tile_names()

    refs = tiles[:4]
    AUTS = tiles[4:]

    # Split the list into XX and YY lists
    refs_XX = refs[::2]
    refs_YY = refs[1::2]

    AUTS_XX = AUTS[::2]
    AUTS_YY = AUTS[1::2]

    # Create a list of pairs of tiles to be aligned.
    # Essentially all possible combinations of Refs, AUTS
    align_pairs = []
    for pair in list(product(refs_XX, AUTS_XX)):
        align_pairs.append(pair)

    for pair in list(product(refs_YY, AUTS_YY)):
        align_pairs.append(pair)

    # Time stuff to help traverse data dir tree.
    t_start = datetime.strptime(start_date, "%Y-%m-%d")
    t_stop = datetime.strptime(stop_date, "%Y-%m-%d")
    n_days = (t_stop - t_start).days

    dates = []
    date_time = []

    for i in range(n_days + 1):
        day = t_start + timedelta(days=i)
        date = day.strftime("%Y-%m-%d")
        dates.append(date)
        d_t = []

        for j in range(48):
            t_delta = datetime.strptime(date, "%Y-%m-%d") + timedelta(minutes=30 * j)
            d_time = t_delta.strftime("%Y-%m-%d-%H:%M")
            d_t.append(d_time)

        date_time.append(d_t)

    data_dir = Path(data_dir)
    out_dir = Path(out_dir)

    for pair in align_pairs:
        for d in range(len(dates)):
            ref = pair[0]
            aut = pair[1]

            ref_path = data_dir / ref / dates[d]
            aut_path = data_dir / aut / dates[d]

            # Parallization magic happens here
            with concurrent.futures.ProcessPoolExecutor() as executor:
                results = executor.map(
                    save_aligned,
                    repeat(ref),
                    repeat(aut),
                    repeat(ref_path),
                    repeat(aut_path),
                    repeat(savgol_window_1),
                    repeat(savgol_window_2),
                    repeat(polyorder),
                    repeat(interp_type),
                    repeat(interp_freq),
                    repeat(out_dir),
                    repeat(dates[d]),
                    date_time[d],
                )

            for result in results:
                print(result)
