"""
MWA Pointings
-------------
Tools to download metadata of the MWA telescope and extract its observational schedule

"""

import numpy as np
import json, wget, time
from pathlib import Path
from astropy.time import Time


def download_meta(start, stop, num_pages, out_dir):
    """Download MWA metadata from `mwatelescope.org <http://mwatelescope.org/>`_

    :param start: start date in :samp:`isot` format :samp:`YYYY-MM-DDTHH:MM:SS` :class:`~str` 
    :param stop: stop date in :samp:`isot` format :samp:`YYYY-MM-DDTHH:MM:SS` :class:`~str`
    :param num_pages: Each page contains 200 observation. Visit `ws.mwatelescope.org/metadata/find <http://ws.mwatelescope.org/metadata/find>`_ to find the total number of pages :class:`~int`
    :param out_dir: Path to output directory where metadata will be saved :class:`~str`

    :returns:
        MWA metadata json files saved to :samp:`out_dir`

    """

    print("Due to download limits, this will take a while")
    wait = 29  # seconds between downloads
    t = wait * num_pages
    m, _ = divmod(t, 60)
    h, m = divmod(m, 60)
    print(f"ETA: Approximately {h:d}H:{m:02d}M")

    # convert isot time to gps
    start_gps = int(Time(start, format="isot").gps)
    stop_gps = int(Time(stop, format="isot").gps)
    
    mwa_meta_dir = Path(out_dir/mwa_pointings)
    mwa_meta_dir.mkdir(parents=True, exist_okay=True)
    for npg in range(num_pages):
        time.sleep(wait)
        cerberus_url = f"http://ws.mwatelescope.org/metadata/find?mintime={start_gps}&maxtime={stop_gps}&extended=1&page={npg+1}&pretty=1"
        print(f"\nDownloading page {npg+1} of metadata")
        wget.download(cerberus_url, f"{mwa_meta_dir}/page_{npg+1:03d}.json")


def clean_meta_json(meta_dir):
    """Organize json files. Clean up quirks in metadata 

    Clear up confusion in Andrew Williams' custom pointing names
    
    :param meta_dir: Path to root of directory with json metadata files :class:`~str`

    :returns:
        A :class:`~tuple` (start_gps, stop_gps, obs_length, pointings)
        
        - start_gps: :class:`~list` of observation start times in :samp:`gps` format
        - stop_gps: :class:`~list` of observation start times in :samp:`gps` format
        - obs_length: :class:`~list` of observation durations in :samp:`seconds`
        - pointings: :class:`~list` of observation :samp:`pointings`

    """

    start_gps = []
    stop_gps = []
    obs_length = []
    pointings = []

    point_0 = ["All_0", "Zenith_Test"]
    point_2 = ["EOR_Point_2", "EOR_Point_2_Delays0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3_Ch100"]
    point_4 = ["EOR_Point_4", "EOR_Point_4_Delays3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0_Ch100"]

    files = Path(f"{meta_dir}/mwa_pointings").glob("page*.json")

    for f in files:
        with open(f) as table:
            data = json.load(table)
            for i in range(len(data)):
                name = data[i][2]
                pointing = data[i][-1]
                start = data[i][0]

                # Many obs are named None, but their title indicates their true pointing
                # These are mostly obs which Andrew Williams schduled for this experiment
                if name in point_0:
                    pointing = 0
                elif name in point_2:
                    pointing = 2
                elif name in point_4:
                    pointing = 4
                else:
                    pass

                ## Filter obs by passes
                # if pointing in [0, 2, 4]:
                if pointing in list(range(197)):

                    # not last obs
                    if i != (len(data) - 1):
                        # Telescope pointing remains constant till it is changed
                        # So, stop time is the start of next observation
                        stop = data[i + 1][0]
                    else:
                        stop = data[i][1]

                    length = stop - start

                    # if length >= 120:

                    start_gps.append(start)
                    stop_gps.append(stop)
                    obs_length.append(length)
                    pointings.append(pointing)

    return (start_gps, stop_gps, obs_length, pointings)


def combine_pointings(start_gps, stop_gps, obs_length, pointings, out_dir):
    """
    Combine successive observations with same pointing and save to file.
        
    :param start_gps: :class:`~list` of observation start times in :samp:`gps` format
    :param stop_gps: :class:`~list` of observation start times in :samp:`gps` format
    :param obs_length: :class:`~list` of observation durations in :samp:`seconds`
    :param pointings: :class:`~list` of observation :samp:`pointings`
    :param out_dir: Path to output directory where cleaned data will be saved :class:`~str`

    :returns:
        :samp:`mwa_pointing.json` saved to :samp:`out_dir`
        
    """

    # if consecutive obs have same pointing, combine them.
    for i in range(len(start_gps)):

        # Not last observation
        if i != len(start_gps) - 1:

            # If consecutive obs have same pointing
            if (stop_gps[i] == start_gps[i + 1]) and (pointings[i] == pointings[i + 1]):
                start_gps[i + 1] = start_gps[i]
                obs_length[i + 1] += obs_length[i]

                pointings[i] = None

    # Apply mask, couldn't think of a better way to implement this
    good_idx = np.where(np.asarray(pointings) != None)[0]
    start_gps = np.asarray(start_gps)[good_idx].tolist()
    stop_gps = np.asarray(stop_gps)[good_idx].tolist()
    obs_length = np.asarray(obs_length)[good_idx].tolist()
    pointings = np.asarray(pointings)[good_idx].tolist()

    # Create dictionary to be saved to json, for plotting
    pointing_list = {}
    pointing_list["grid_pt"] = pointings
    pointing_list["start_gps"] = start_gps
    pointing_list["stop_gps"] = stop_gps
    pointing_list["obs_length"] = obs_length
    
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    with open(f"{out_dir}/mwa_pointing.json", "w") as outfile:
        json.dump(pointing_list, outfile, indent=4)


if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(
        description="""
            Download pointing metadata
            """
    )

    parser.add_argument(
        "--start",
        metavar="\b",
        help="Start date in isot format. Ex: 2019-10-01T00:00:00",
    )
    parser.add_argument(
        "--stop",
        metavar="\b",
        help="Stop date in isot format. Ex: 2019-10-10T00:00:00",
    )
    parser.add_argument(
        "--num_pgs",
        type=int,
        metavar="\b",
        help="Total number of pages between start and stop, determined from http://ws.mwatelescope.org/metadata/find",
    )
    parser.add_argument(
        "--out_dir", metavar="\b", help="Output directory for json metadata",
    )

    args = parser.parse_args()
    start = args.start
    stop = args.stop
    num_pgs = args.num_pgs
    out_dir = Path(args.out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)

    download_meta(start, stop, num_pgs, out_dir)
