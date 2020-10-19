"""
MWA Pointings
-------------
Tools to download metadata of the MWA telescope and extract its observational schedule

"""
import json
import time
from datetime import datetime, timedelta
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pytz
import seaborn as sns
import wget
from astropy.time import Time
from embers.rf_tools.rf_data import tile_names, time_tree
from matplotlib import pyplot as plt

mpl.use("Agg")


def download_meta(start, stop, num_pages, out_dir, wait):
    """Download MWA metadata from `mwatelescope.org <http://mwatelescope.org/>`_

    :param start: start date in :samp:`isot` format :samp:`YYYY-MM-DDTHH:MM:SS` :class:`~str`
    :param stop: stop date in :samp:`isot` format :samp:`YYYY-MM-DDTHH:MM:SS` :class:`~str`
    :param num_pages: Each page contains 200 observation. Visit `ws.mwatelescope.org/metadata/find <http://ws.mwatelescope.org/metadata/find>`_ to find the total number of pages :class:`~int`
    :param out_dir: Path to output directory where metadata will be saved :class:`~str`
    :param wait: Time to sleep between downloads so as not to overload servers. Default=29

    :returns:
        MWA metadata json files saved to :samp:`out_dir`

    """

    print("Downloading MWA metadata")
    print("Due to download limits, this will take a while")
    t = wait * num_pages
    m, _ = divmod(t, 60)
    h, m = divmod(m, 60)
    print(f"ETA: Approximately {h:d}H:{m:02d}M")

    # convert isot time to gps
    start_gps = int(Time(start, format="isot").gps)
    stop_gps = int(Time(stop, format="isot").gps)

    mwa_meta_dir = Path(f"{out_dir}/mwa_pointings")
    mwa_meta_dir.mkdir(parents=True, exist_ok=True)
    for npg in range(num_pages):
        time.sleep(wait)
        cerberus_url = f"http://ws.mwatelescope.org/metadata/find?mintime={start_gps}&maxtime={stop_gps}&extended=1&page={npg+1}&pretty=1"
        print(f"\nDownloading page {npg+1}/{num_pages} of metadata")
        wget.download(cerberus_url, f"{mwa_meta_dir}/page_{npg+1:03d}.json")


def clean_meta_json(out_dir):
    """Organize json files. Clean up quirks in metadata

    Clear up confusion in Andrew Williams' custom pointing names

    :param out_dir: Path to root of directory with json metadata files :class:`~str`

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

    files = Path(f"{out_dir}/mwa_pointings").glob("page*.json")

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

                # Filter obs by passes
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
    with open(f"{out_dir}/mwa_pointings.json", "w") as outfile:
        json.dump(pointing_list, outfile, indent=4)


def point_integration(out_dir):
    """Calculate total integration at each pointing

    :param out_dir: Path to directory where :samp:`mwa_pointings.json` is saved

    :returns:
        A :class:`~tuple`

        - pointings : :class:`~list` of MWA pointings
        - int_hours: :class:`~list` of total integration, at each pointing, in hours

    """

    with open(f"{out_dir}/mwa_pointings.json") as data:
        pointings = json.load(data)
        grid_pt = pointings["grid_pt"]
        obs_length = pointings["obs_length"]

    # find unique pointings.
    unique_pointings = list(set(grid_pt))
    unique_pointings = sorted(unique_pointings)

    pointings = []
    integrations = []

    for i in unique_pointings:
        pointings.append(i)
        time = 0
        for j in range(len(grid_pt)):
            if grid_pt[j] == i:
                time += obs_length[j]
        integrations.append(time)

    int_hours = np.asarray(integrations) / (60 * 60)

    return (pointings, int_hours)


def pointing_hist(pointings, int_hours, time_thresh, out_dir):
    """
    Plot a histogram of pointing integration

    Many pointings can have very low total integration, use the :samp:`time_thresh` argument to exclude low integration pointings.

    :param pointings: :class:`~list` of MWA pointings
    :param int_hours: :class:`~list` of total integration, at each pointing, in hours
    :param time_thresh: minimum integration time to be included in histogram, in hours :class:`~int`
    :param out_dir: Path to directory where :samp:`mwa_pointings.json` is saved

    :returns:
        Pointing histogram plot saved to :samp:`out_dir`

    """

    # Filter out pointings below time_thresh
    time_point = []
    point = []

    for i in range(len(int_hours)):
        if int_hours[i] >= time_thresh:
            point.append(pointings[i])
            time_point.append(int_hours[i])

    x = range(len(time_point))
    leg = [int(i) for i in time_point]

    plt.style.use("seaborn")
    _, ax = plt.subplots(figsize=(8, 6))
    pal = sns.cubehelix_palette(
        len(time_point), start=0.4, rot=-0.5, dark=0.4, reverse=True
    )
    barplot = plt.bar(x, time_point, color=sns.color_palette(pal))

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

    plt.xticks(x, point)
    plt.ylabel("Hours")
    plt.xlim(-0.7, len(time_point) - 0.3)
    plt.xlabel("MWA Grid Pointing Number")
    plt.title("Integration at MWA Grid Pointings")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/pointing_integration.png")
    print(f"Pointing integration plot saved to {out_dir}")


def rf_obs_times(start_date, stop_date, time_zone):
    """Generate start & end times of 30 minuts rf observations in local, unix, gps formats

    :param start_date: in :samp:`YYYY-MM-DD` format :class:`~str`
    :param stop_date: in :samp:`YYYY-MM-DD` format :class:`~str`
    :param time_zone: A :class:`~str` representing a :samp:`pytz` `timezones <https://gist.github.com/heyalexej/8bf688fd67d7199be4a1682b3eec7568>`_.

    :returns:
        A :class:`~tuple` (obs_time, obs_unix, obs_unix_end, obs_gps, obs_gps_end)

        - obs_time: :class:`~list` of start times of 30 min obs in :samp:`YYYY-MM-DD-HH:MM` format
        - obs_gps: :class:`~list` of start times of 30 min obs in :samp:`gps` format
        - obs_gps_end: :class:`~list` of end times of 30 min obs in :samp:`gps` format

    """

    # The time input is in local time. As in Austraila/Perth
    local = pytz.timezone(time_zone)

    t_start = datetime.strptime(start_date, "%Y-%m-%d")
    t_stop = datetime.strptime(stop_date, "%Y-%m-%d")

    # Number of days that the date range spans
    n_days = (t_stop - t_start).days

    # YYYY-MM-DD-HH:MM format
    obs_time = []

    # Start of half hour obs in unix time
    obs_unix = []

    # End of half hour obs in unix time
    obs_unix_end = []

    # The +1 makes the date ranges inclusive
    for i in range(n_days + 1):
        day = t_start + timedelta(days=i)
        date = day.strftime("%Y-%m-%d")

        # loop over 48x30 minute obs in a day
        for j in range(48):
            t_delta = datetime.strptime(date, "%Y-%m-%d") + timedelta(minutes=30 * j)
            # convert t_delta to a readable string YYYY-MM-DD-HH:MM
            d_time = t_delta.strftime("%Y-%m-%d-%H:%M")

            # Convert from naive local time to utc aware time
            utc_delta = local.localize(t_delta, is_dst=None).astimezone(pytz.utc)

            # convert to a unix timestamp, used within rf explorer data files
            utc_unix = utc_delta.timestamp()
            # time at end of half hour window
            utc_unix_end = utc_unix + (30 * 60)

            obs_time.append(d_time)
            obs_unix.append(utc_unix)
            obs_unix_end.append(utc_unix_end)

    # Start and end of 30 min obs in gps time
    # Round to nearest int
    obs_gps = np.rint(Time(obs_unix, format="unix").gps)
    obs_gps_end = np.rint(Time(obs_unix_end, format="unix").gps)

    return (obs_time, obs_gps, obs_gps_end)


def obs_pointings(start, stop, time_zone, out_dir):
    """
    Classify the pointing of each :samp:`rf_obs`

    Loop over all rf observations within a date interval and determine whether
    the 30 minute period had more that a 60% majority at a single pointing. If
    it does, the rf observation is saved to an appropriate list. Save the pointing
    data to :samp:`obs_pointings.json` in the :samp:`out_dir`.

    :param start: in :samp:`YYYY-MM-DD` format :class:`~str`
    :param stop: in :samp:`YYYY-MM-DD` format :class:`~str`
    :param time_zone: A :class:`~str` representing a :samp:`pytz` `timezones <https://gist.github.com/heyalexej/8bf688fd67d7199be4a1682b3eec7568>`_.
    :param out_dir: Path to directory where :samp:`mwa_pointings.json` is saved

    :returns:
        :samp:`obs_pointings.json` saved to :samp:`out_dir`

    """

    point_0 = []
    point_2 = []
    point_4 = []
    point_41 = []

    obs_time, obs_gps, obs_gps_end = rf_obs_times(start, stop, time_zone)

    with open(f"{out_dir}/mwa_pointings.json") as table:
        data = json.load(table)
        pointings = data["grid_pt"]
        start_gps = data["start_gps"]
        stop_gps = data["stop_gps"]

    # loop over each 30 min rf_obs
    for i in range(len(obs_time)):

        # loop over mwa pointing observations from
        # ultimate_pointing_times.json
        for j in range(len(start_gps)):

            # start and stop gps times of MWA pointing observations
            meta_start = start_gps[j]
            meta_stop = stop_gps[j]

            # I. MWA obs starts before rf_obs, stops within rf_obs
            if (
                meta_start < obs_gps[i]
                and meta_stop > obs_gps[i]
                and meta_stop <= obs_gps_end[i]
            ):
                occu = (meta_stop - obs_gps[i]) / 1800

            # II. MWA obs completely within rf_obs
            elif meta_start >= obs_gps[i] and meta_stop <= obs_gps_end[i]:
                occu = (meta_stop - meta_start) / 1800

            # III. MWA obs starts within rf obs and ends after
            elif (
                meta_start >= obs_gps[i]
                and meta_start < obs_gps_end[i]
                and meta_stop > obs_gps_end[i]
            ):
                occu = (obs_gps_end[i] - meta_start) / 1800

            # IV. MWA obs starts before and ends after rf obs
            elif meta_start < obs_gps[i] and meta_stop > obs_gps_end[i]:
                occu = (obs_gps_end[i] - obs_gps[i]) / 1800

            # V. No common time between MWA obs and rf obs
            else:
                occu = 0

            if occu != 0 and occu >= 0.6:
                if pointings[j] == 0:
                    point_0.append(obs_time[i])
                elif pointings[j] == 2:
                    point_2.append(obs_time[i])
                elif pointings[j] == 4:
                    point_4.append(obs_time[i])
                elif pointings[j] == 41:
                    point_41.append(obs_time[i])
                # should there be an else statement here?

    # Create dictionary to be saved to json
    obs_pointings = {}

    obs_pointings["start_date"] = start
    obs_pointings["stop_date"] = stop
    obs_pointings["point_0"] = point_0
    obs_pointings["point_2"] = point_2
    obs_pointings["point_4"] = point_4
    obs_pointings["point_41"] = point_41

    with open(f"{out_dir}/obs_pointings.json", "w") as outfile:
        json.dump(obs_pointings, outfile, indent=4)


def tile_integration(out_dir, rf_dir):
    """Calculate total integration at multiple pointings for all tiles

    :param out_dir: Path to root of directory where mwa metadata is saved :class:`~str`
    :param rf_dir: Path to root of directory with rf data files :class:`~str`

    :returns:
        A :class:`~dict` with keys being :samp:`tile_names` and values being :class:`~list` of integration hours at pointings 0, 2, 4 & 41

    """

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

    tile_ints = {}
    for tile in tiles:

        # increment this by 0.5, for every obs which is in point_list.
        # 0.5 = 30 min
        p_0 = 0
        p_2 = 0
        p_4 = 0
        p_41 = 0

        # dates: list of days
        # date_time = list of 30 min observation windows
        dates, time_stamps = time_tree(start_date, stop_date)
        for day in range(len(dates)):
            for ts in range(len(time_stamps[day])):
                f_name = f"{tile}_{time_stamps[day][ts]}.txt"
                f_path = Path(f"{rf_dir}/{tile}/{dates[day]}/{f_name}")
                if f_path.is_file():
                    if time_stamps[day][ts] in point_0:
                        p_0 += 0.5
                    elif time_stamps[day][ts] in point_2:
                        p_2 += 0.5
                    elif time_stamps[day][ts] in point_4:
                        p_4 += 0.5
                    elif time_stamps[day][ts] in point_41:
                        p_41 += 0.5
                    # should there be an else statement here?

        tile_ints[f"{tile}"] = [p_0, p_2, p_4, p_41]

    return tile_ints


def plt_hist_array(tile_ints, out_dir):
    """A massive grid of histograms with a subplot for pointing integration of each tile.

    :param tile_ints: :class:`~dict` from :func:`~embers.mwa_utils.mwa_pointings.tile_integration`
    :param out_dir: Path to output directory :class:`~str`

    :returns:
        :samp:`tiles_pointing_integration.png` saved to :samp:`out_dir`

    """

    plt.style.use("seaborn")
    fig = plt.figure(figsize=(14, 10))

    # Tile names
    tiles = tile_names()

    def tile_pointing_hist(a, b, c, fig, int_dict=None, tile=None):

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

    tile_pointing_hist(4, 8, 1, fig, int_dict=tile_ints, tile=tiles[0])
    tile_pointing_hist(4, 8, 2, fig, int_dict=tile_ints, tile=tiles[4])
    tile_pointing_hist(4, 8, 3, fig, int_dict=tile_ints, tile=tiles[5])
    tile_pointing_hist(4, 8, 4, fig, int_dict=tile_ints, tile=tiles[6])
    tile_pointing_hist(4, 8, 5, fig, int_dict=tile_ints, tile=tiles[7])
    tile_pointing_hist(4, 8, 6, fig, int_dict=tile_ints, tile=tiles[8])
    tile_pointing_hist(4, 8, 7, fig, int_dict=tile_ints, tile=tiles[9])
    tile_pointing_hist(4, 8, 8, fig, int_dict=tile_ints, tile=tiles[10])
    tile_pointing_hist(4, 8, 9, fig, int_dict=tile_ints, tile=tiles[1])
    tile_pointing_hist(4, 8, 10, fig, int_dict=tile_ints, tile=tiles[11])
    tile_pointing_hist(4, 8, 11, fig, int_dict=tile_ints, tile=tiles[12])
    tile_pointing_hist(4, 8, 12, fig, int_dict=tile_ints, tile=tiles[13])
    tile_pointing_hist(4, 8, 13, fig, int_dict=tile_ints, tile=tiles[14])
    tile_pointing_hist(4, 8, 14, fig, int_dict=tile_ints, tile=tiles[15])
    tile_pointing_hist(4, 8, 15, fig, int_dict=tile_ints, tile=tiles[16])
    tile_pointing_hist(4, 8, 16, fig, int_dict=tile_ints, tile=tiles[17])
    tile_pointing_hist(4, 8, 17, fig, int_dict=tile_ints, tile=tiles[2])
    tile_pointing_hist(4, 8, 18, fig, int_dict=tile_ints, tile=tiles[18])
    tile_pointing_hist(4, 8, 19, fig, int_dict=tile_ints, tile=tiles[19])
    tile_pointing_hist(4, 8, 20, fig, int_dict=tile_ints, tile=tiles[20])
    tile_pointing_hist(4, 8, 21, fig, int_dict=tile_ints, tile=tiles[21])
    tile_pointing_hist(4, 8, 22, fig, int_dict=tile_ints, tile=tiles[22])
    tile_pointing_hist(4, 8, 23, fig, int_dict=tile_ints, tile=tiles[23])
    tile_pointing_hist(4, 8, 24, fig, int_dict=tile_ints, tile=tiles[24])
    tile_pointing_hist(4, 8, 25, fig, int_dict=tile_ints, tile=tiles[3])
    tile_pointing_hist(4, 8, 26, fig, int_dict=tile_ints, tile=tiles[25])
    tile_pointing_hist(4, 8, 27, fig, int_dict=tile_ints, tile=tiles[26])
    tile_pointing_hist(4, 8, 28, fig, int_dict=tile_ints, tile=tiles[27])
    tile_pointing_hist(4, 8, 29, fig, int_dict=tile_ints, tile=tiles[28])
    tile_pointing_hist(4, 8, 30, fig, int_dict=tile_ints, tile=tiles[29])
    tile_pointing_hist(4, 8, 31, fig, int_dict=tile_ints, tile=tiles[30])
    tile_pointing_hist(4, 8, 32, fig, int_dict=tile_ints, tile=tiles[31])

    plt.tight_layout()
    fig.savefig(f"{out_dir}/tiles_pointing_integration.png")
    print(f"Tile integration plot saved to {out_dir}")


def mwa_point_meta(
    start, stop, num_pages, time_thresh, time_zone, rf_dir, out_dir, wait=29
):
    """
    Download mwa pointing metadata, sort and parse it, and create diagonistic plots,

    :param start: start date in :samp:`isot` format :samp:`YYYY-MM-DDTHH:MM:SS` :class:`~str`
    :param stop: stop date in :samp:`isot` format :samp:`YYYY-MM-DDTHH:MM:SS` :class:`~str`
    :param num_pages: Each page contains 200 observation. Visit `ws.mwatelescope.org/metadata/find <http://ws.mwatelescope.org/metadata/find>`_ to find the total number of pages :class:`~int`
    :param time_thresh: minimum integration time to be included in histogram, in hours :class:`~int`
    :param time_zone: A :class:`~str` representing :samp:`pytz` `timezones <https://gist.github.com/heyalexej/8bf688fd67d7199be4a1682b3eec7568>`_.
    :param rf_dir: Path to root of directory with rf data files :class:`~str`
    :param out_dir: Path to output directory where metadata will be saved :class:`~str`
    :param wait: Time to sleep between downloads so as not to overload servers. Default=29

    :returns:
        Data products saved to :samp:`out_dir`

    """

    # Download pointing metadata
    download_meta(start, stop, num_pages, out_dir, wait)
    print("\nMetadata download complete")

    # Organize and combine metadata
    start_gps, stop_gps, obs_length, pointings = clean_meta_json(out_dir)
    combine_pointings(start_gps, stop_gps, obs_length, pointings, out_dir)

    # Compute pointing integration and plot histogram
    pointings, int_hours = point_integration(out_dir)
    pointing_hist(pointings, int_hours, time_thresh, out_dir)

    # Determine pointing of each 30 minute observation
    obs_pointings(
        start, stop, time_zone, out_dir,
    )

    # Pointing integration for each tile
    tile_ints = tile_integration(out_dir, rf_dir)
    plt_hist_array(tile_ints, out_dir)
