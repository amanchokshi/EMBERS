"""
Satellite Channels
------------------

A set of tools to determine the transmission channel of various satellites in 30 minute
observation by using rf data in conjunctin with chronological satellite ephemeris data.

"""

import re
import json
import argparse
import numpy as np
import seaborn as sns
from pathlib import Path
import concurrent.futures
from itertools import repeat
import matplotlib.pylab as pl
import matplotlib.pyplot as plt
from scipy.stats import median_absolute_deviation as mad

from embers.rf_tools.colormaps import spectral
from embers.rf_tools.rf_data import time_tree


def read_ref_aligned(ali_file=None):
    """Read aligned reference data from :func:`~embers.rf_tools.align_data.save_aligned` :samp:`npz` file
    
    :param ali_file: path to a :func:`~embers.rf_tools.align_data.save_aligned` :samp:`npz` file :class:`~str`

    :returns:
        A :class:`~tuple` (power, times)

        - power: a smoothed reference power array :class:`~numpy.ndarry`
        - times: a regular time array :class:`~numpy.ndarry`

    """

    paired_data = np.load(ali_file, allow_pickle=True)

    power = paired_data["ref_ali"]
    times = paired_data["time_array"]

    return (power, times)


def noise_floor(sat_thresh, noi_thresh, power):
    """Computes the noise floor of a rf power array
    
    Exclude channels with signal above :samp:`sat_thresh` multiplied by :samp:`standard deviation` of power array.
    The Median Absolute Deviation :samp:`MAD` is used to quantify the noise level of the remaining
    channels. The noise floor :samp:`noi_thresh` is defined to be the :samp:`median` of noisy data + :samp:`noi_thresh` multiplied by the
    :samp:`MAD` of noisy data.

    :param sat_thresh: An integer multiple of standard deviation of rf power array, used to exclude channels with potential satellites. :class:`~int` 
    :param noi_thresh: An integer multiple of the noisy data MAD, used to compute a noise floor. :class:`~int`
    :param power: Rf power array :class:`~numpy.ndarry`

    :returns:
        noise_threshold: The power level of the noise floor in dBm :class:`~int`

    """

    # compute the standard deviation of data, and use it to identify occupied channels
    σ = np.std(power)

    # Any channel with a max power >= σ has a satellite
    sat_cut = sat_thresh * σ
    chans_pow_max = np.amax(power, axis=0)

    # Exclude the channels with sats, to only have noise data
    noise_chans = np.where(chans_pow_max < sat_cut)[0]
    noise_data = power[:, noise_chans]

    # noise median, noise mad, noise threshold = μ + 3*σ
    μ_noise = np.median(noise_data)
    σ_noise = mad(noise_data, axis=None)
    noise_threshold = μ_noise + noi_thresh * σ_noise

    return noise_threshold


def time_filter(s_rise, s_set, times):
    """Determine indices of time array when a satellite is above the horizon.

    Isolate the portion of a rf power array where a satellite is above the
    horizon using the rise and set times of a satellite's ephemeris. This 
    function returns a pair of indices which can be used to slice the rf 
    power and times arrays to precisely only include the satellite.

    :param s_rise: satellite rise time from ephemeris :class:`~float`
    :param s_set: satellite set time from ephemeris :class:`~float`
    :param times: time array corresponding to the rf power array :class:`~numpy.ndarry`

    :returns:
        intvl: :samp:`None` if satellite is not above the horizon within the :samp:`times` array
        intvl: [i_0, i_1], pair of indices of :samp:`times` array, when satellite is above the horizon
    
    """

    # I. sat rises before times, sets within times window
    if s_rise < times[0] and s_set > times[0] and s_set <= times[-1]:
        i_0 = np.where(times == times[0])[0][0]
        i_1 = np.where(times == s_set)[0][0]
        intvl = [i_0, i_1]

    # II. sat rises and sets within times
    elif s_rise >= times[0] and s_set <= times[-1]:
        i_0 = np.where(times == s_rise)[0][0]
        i_1 = np.where(times == s_set)[0][0]
        intvl = [i_0, i_1]

    # III. sat rises within times, and sets after
    elif s_rise >= times[0] and s_rise < times[-1] and s_set > times[-1]:
        i_0 = np.where(times == s_rise)[0][0]
        i_1 = np.where(times == times[-1])[0][0]
        intvl = [i_0, i_1]

    # IV. sat rises before times and sets after
    elif s_rise < times[0] and s_set > times[-1]:
        i_0 = np.where(times == times[0])[0][0]
        i_1 = np.where(times == times[-1])[0][0]
        intvl = [i_0, i_1]

    # V. sat completely out of times. Could be on either side
    else:
        intvl = None

    # intvl = interval
    return intvl


def plt_window_chans(power, sat_id, start, stop, cmap, chs=None, good_ch=None):
    """Waterfall plot with sat window and occupied channels highlighted.
    
    :param power: Rf power array :class:`~numpy.ndarry`
    :param sat_id: Norad catalogue ID :class:`~str`
    :param start: Index of power array when :samp:`sat_id` is above the horizon :class:`~int`
    :param stop: Index of power array when :samp:`sat_id` is above the horizon :class:`~int`
    :param cmap: Colormap for plotting waterfall :class:`~matplotlib.colors.ListedColormap`
    :param chs: Occupied channels :class:`~list` of :class:`~int`
    :param good_ch: Most probable channel :class:`~int`
    
    :returns:
        - plt - :func:`~matplotlib.pyplot.plot` object

    """

    plt.rcParams.update(plt.rcParamsDefault)
    plt.style.use("dark_background")
    fig = plt.figure(figsize=(7, 10))
    ax = fig.add_axes([0.12, 0.1, 0.72, 0.85])
    im = ax.imshow(power, interpolation="none", cmap=cmap)
    cax = fig.add_axes([0.88, 0.1, 0.03, 0.85])
    fig.colorbar(im, cax=cax)
    ax.set_aspect("auto")
    ax.set_title(f"Waterfall Plot: {sat_id} in {chs}")

    ax.set_xlabel("Freq Channel")
    ax.set_ylabel("Time Step")

    # horizontal highlight of satellite time window
    ax.axhspan(start, stop, alpha=0.1, color="white")

    if chs != None:
        # vertical highlight of possible channels
        for ch in chs:
            ax.axvspan(ch - 0.5, ch + 0.4, alpha=0.2, color="white")

    if good_ch != None:
        # highlight good channel
        ax.axvspan(
            good_ch - 0.5,
            good_ch + 0.4,
            alpha=0.6,
            edgecolor="#a7ff83",
            facecolor=None,
            fill=False,
        )

    return plt


def plt_channel(
    times, channel_power, pow_med, chan_num, y_range, noi_thresh, pow_thresh
):
    """Plot power in channel, with various thresholds
    
    :param times: times 1D array :class:`~numpy.ndarry`
    :param channel_power: Rf channel power 1D array :class:`~numpy.ndarry`
    :param pow_med: Median of rf power array :class:`~float`
    :param chan_num: Channel number :class:`~str`
    :param y_range: Plot min, max yrange list :class:`~list`
    :param noi_thresh: Noise floor in :samp:`dBm`. :class:`~float`
    :param pow_thresh: Power threshold in :samp:`dBm` :class:`~float` 
        
    :returns:
        - plt - :func:`~matplotlib.pyplot.plot` object

    """
    plt.rcParams.update(plt.rcParamsDefault)
    plt.style.use("seaborn")

    # plt channel power
    plt.plot(
        times,
        channel_power,
        linestyle="-",
        linewidth=2,
        alpha=1.0,
        color="#db3751",
        label="Data",
    )
    plt.fill_between(
        times,
        np.repeat(y_range[0], len(channel_power)),
        channel_power,
        color="#db3751",
        alpha=0.7,
    )

    plt.axhline(
        pow_thresh,
        alpha=1.0,
        linestyle="-",
        linewidth=2,
        color="#fba95f",
        label=f"Power Cut: {pow_thresh:.2f} dBm",
    )

    plt.axhline(
        noi_thresh,
        linestyle="-",
        linewidth=2,
        color="#5cb7a9",
        label=f"Noise Cut: {noi_thresh:.2f} dBm",
    )

    plt.ylim(y_range)
    plt.xlim([times[0], times[-1]])
    plt.ylabel("Power [dBm]")
    plt.xlabel("Time [s]")
    plt.title(f"Satellite Pass in Channel: [{chan_num}]")
    plt.tight_layout()
    leg = plt.legend(frameon=True)
    leg.get_frame().set_facecolor("grey")
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)

    return plt


def plt_sats(ids, alt, az, timestamp):
    """Polar plot of satellite passes in a 30 minute observation
    
    :param ids: :class:`~list` of Norad catalogue IDs
    :param alt: :class:`~list` of :samp:`Altitudes` values
    :param az: :class:`~list` of :samp:`~Azimuth` values
    :param timestamp: Time at start of 30 minute observation in :samp:`YYYY-MM-DD-HH:MM` format :class:`~str`
    
    :returns:
        - plt - :func:`~matplotlib.pyplot.plot` object

    """

    plt.rcParams.update(plt.rcParamsDefault)
    plt.style.use("seaborn")
    figure = plt.figure(figsize=(7, 6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0, 30, 60, 90], angle=22)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_title(f"Satellite passes in {timestamp}", y=1.05)
    ax.grid(color="grey", linewidth=1.6, alpha=0.5)

    colors = pl.cm.Spectral(np.linspace(0.17, 0.9, len(alt)))

    for i in range(len(alt)):
        plt.plot(
            az[i],
            alt[i],
            "-",
            linewidth=2.4,
            alpha=0.8,
            color=colors[i],
            label=f"{ids[i]}",
        )
        plt.legend()

    leg = plt.legend(
        frameon=True,
        bbox_to_anchor=(0.28, 1.0, 1.0, -0.95),
        loc="center right",
        title="NoradID",
    )
    leg.get_frame().set_facecolor("grey")
    leg.get_frame().set_alpha(0.4)
    for l in leg.legendHandles:
        l.set_alpha(1)

    plt.tight_layout()

    return plt


def good_chans(
    ali_file,
    chrono_file,
    sat_id,
    sat_thresh,
    noi_thresh,
    pow_thresh,
    occ_thresh,
    timestamp,
    out_dir,
    plots=None,
):
    """Determine the channels a satellite could occupy, in a 30 minute observation
    
    Ephemeris from :samp:`chrono_file` is used to select a temoral :samp:`window`
    of the rf power array, within which the satellite is above the horizon. Looping
    through the frequency channels, a :samp:`noi_thresh`, :samp:`pow_thresh`, :samp:`occ_thresh`
    are used to identify possible channels occupied by the :samp:`sat_id`. If more than one
    channel passes the three thresholds, the channel with the highest window occupancy is 
    selected.

    .. code-block:: python
        
        from embers.sat_utils.sat_channels import good_chans
        
        ali_file = "~/embers_out/rf0XX_S06XX_2019-10-10-02:30_aligned.npz"
        chrono_file = "~/embers_out/2019-10-10-02:30.json"
        sat_id = "44387"
        sat_thresh = 1
        noi_thresh = 3
        pow_thresh = 20
        occ_thresh = 0.80
        timestamp = "2019-10-10-02:30"
        out_dir = "./embers_out"
        plots = True

        good_chan = good_chans(
                        ali_file,
                        chrono_file,
                        sat_id,
                        sat_thresh,
                        noi_thresh,
                        pow_thresh,
                        occ_thresh,
                        timestamp,
                        out_dir,
                        plots=plots)

        print(good_chan)
        >>> 59

    :param ali_file: Path to a :samp:`npz` aligned file from :func:`~embers.rf_data.align_data.save_aligned` :class:`~str`
    :param chrono_file: Path to chrono ephem json file from :func:`~embers.sat_utils.chrono_ephem.save_chrono_ephem` :class:`~str`
    :param sat_id: Norad catalogue ID :class:`~str`
    :param sat_thresh: Satellite threshold from :func:`~embers.sat_utils.sat_channels.noise_threshold` :class:`~int`
    :param noi_thresh: Noise threshold from :func:`~embers.sat_utils.sat_channels.noise_threshold` :class:`~int`
    :param occ_thresh: Window occupation threshold. Minimum fractional signal above the noise floor in window :class:`~float`
    :param timestamp: Time at start of observation in format :samp:`YYYY-MM-DD-HH:MM` :class:`~str`
    :param out_dir: Path to output directory to save plots :class:`~str`
    :param plots: If :samp:`True`, disagnostic plots are generated and saved to :samp:`out_dir`

    :returns:
        - good_chan: The channel number of most probable channel for :samp:`sat_id`
        - good_chan may be :samp:`None`, if no possible channels were identified

    """

    spec, _ = spectral()

    power, times = read_ref_aligned(ali_file=ali_file)
    p_med = np.median(power)

    # Determine noise threshold
    noise_threshold = noise_floor(sat_thresh, noi_thresh, power)

    with open(chrono_file) as chrono:
        chrono_ephem = json.load(chrono)

        # All satellites in chrono ephem json file
        norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(len(chrono_ephem))]

        # Index of sat_id in chrono ephem json file
        norad_index = norad_list.index(sat_id)

        # Extract ephemeris for sat_id
        norad_ephem = chrono_ephem[norad_index]

        rise_ephem = norad_ephem["time_array"][0]
        set_ephem = norad_ephem["time_array"][-1]

        # indices of times, when sat rose and set
        intvl = time_filter(rise_ephem, set_ephem, np.asarray(times))

        # If sat was in the sky in the 30 min obs
        if intvl != None:

            # Window start, stop
            w_start, w_stop = intvl
            window_len = w_stop - w_start + 1

            # Slice the power/times arrays to the times of sat pass
            power_c = power[w_start : w_stop + 1, :]
            times_c = times[w_start : w_stop + 1]

            # possible identified channels
            possible_chans = []
            # window occupancy of possible channels
            occu_list = []

            # Loop over every channel
            for s_chan in range(len(power_c[0])):

                channel_power = power_c[:, s_chan]

                # Percentage of signal occupancy above noise threshold
                window_occupancy = (np.where(channel_power >= noise_threshold))[
                    0
                ].size / window_len

                # Power threshold below which satellites aren't counted
                # Only continue if there is signal for more than 80% of satellite pass
                if (
                    max(channel_power) >= p_med + pow_thresh
                    and occ_thresh <= window_occupancy < 1.00
                ):

                    # If satellite window begins from the start of the observation
                    if times[0] == times_c[0]:
                        # last 10 data points (10 seconds) are below the noise threshold
                        if (
                            all(p < noise_threshold for p in channel_power[-11:-1])
                        ) is True:
                            occu_list.append(window_occupancy)
                            possible_chans.append(s_chan)

                            if plots is True:
                                plt = plt_channel(
                                    times_c,
                                    channel_power,
                                    p_med,
                                    s_chan,
                                    [
                                        np.amin(channel_power) - 1,
                                        np.amax(channel_power) + 1,
                                    ],
                                    noise_threshold,
                                    pow_thresh + p_med,
                                )
                                date = re.search(r"\d{4}.\d{2}.\d{2}", timestamp)[0]
                                plt_dir = Path(f"{out_dir}/plots/{date}/{timestamp}")
                                plt_dir.mkdir(parents=True, exist_ok=True)
                                plt.savefig(
                                    f"{plt_dir}/{sat_id}_channel_{s_chan}_{window_occupancy}.png"
                                )
                                plt.close()

                    # if the satellite window ends and the end of the observation
                    elif times[-1] == times_c[-1]:

                        # first 10 data points (10 seconds) are below the noise threshold
                        if (
                            all(p < noise_threshold for p in channel_power[:10])
                        ) is True:
                            occu_list.append(window_occupancy)
                            possible_chans.append(s_chan)

                            if plots is True:
                                plt = plt_channel(
                                    times_c,
                                    channel_power,
                                    p_med,
                                    s_chan,
                                    [
                                        np.amin(channel_power) - 1,
                                        np.amax(channel_power) + 1,
                                    ],
                                    noise_threshold,
                                    pow_thresh + p_med,
                                )
                                date = re.search(r"\d{4}.\d{2}.\d{2}", timestamp)[0]
                                plt_dir = Path(f"{out_dir}/plots/{date}/{timestamp}")
                                plt_dir.mkdir(parents=True, exist_ok=True)
                                plt.savefig(
                                    f"{plt_dir}/{sat_id}_channel_{s_chan}_{window_occupancy}.png"
                                )
                                plt.close()

                    # satellite window completely within the observation
                    else:
                        # first and last 10 data points (10 seconds) are below the noise threshold
                        if (
                            all(p < noise_threshold for p in channel_power[:10])
                            and all(p < noise_threshold for p in channel_power[-11:-1])
                        ) is True:
                            occu_list.append(window_occupancy)
                            possible_chans.append(s_chan)

                            if plots is True:
                                plt = plt_channel(
                                    times_c,
                                    channel_power,
                                    p_med,
                                    s_chan,
                                    [
                                        np.amin(channel_power) - 1,
                                        np.amax(channel_power) + 1,
                                    ],
                                    noise_threshold,
                                    pow_thresh + p_med,
                                )
                                date = re.search(r"\d{4}.\d{2}.\d{2}", timestamp)[0]
                                plt_dir = Path(f"{out_dir}/plots/{date}/{timestamp}")
                                plt_dir.mkdir(parents=True, exist_ok=True)
                                plt.savefig(
                                    f"{plt_dir}/{sat_id}_channel_{s_chan}_{window_occupancy}.png"
                                )
                                plt.close()

            # If channels are identified in the 30 min obs
            n_chans = len(possible_chans)
            if n_chans > 0:

                # The most probable channel is one with the highest occupation
                good_chan = possible_chans[occu_list.index(max(occu_list))]

                if plots is True:
                    plt = plt_window_chans(
                        power,
                        sat_id,
                        w_start,
                        w_stop,
                        spec,
                        chs=possible_chans,
                        good_ch=good_chan,
                    )
                    date = re.search(r"\d{4}.\d{2}.\d{2}", timestamp)[0]
                    plt_dir = Path(f"{out_dir}/plots/{date}/{timestamp}")
                    plt_dir.mkdir(parents=True, exist_ok=True)
                    plt.savefig(f"{plt_dir}/{sat_id}_waterfall_{good_chan}.png")
                    plt.close()

                return good_chan

            else:
                if plots == "True":
                    plt = plt_window_chans(power, sat_id, w_start, w_stop, spec)
                    date = re.search(r"\d{4}.\d{2}.\d{2}", timestamp)[0]
                    plt_dir = Path(f"{out_dir}/plots/{date}/{timestamp}")
                    plt_dir.mkdir(parents=True, exist_ok=True)
                    plt.savefig(f"{plt_dir}/{sat_id}_waterfall_window.png")
                    plt.close()

                return None

        else:
            return None


def window_chan_map(
    ali_dir,
    chrono_dir,
    sat_thresh,
    noi_thresh,
    pow_thresh,
    occ_thresh,
    out_dir,
    plots,
    plt_dir,
    obs_stamp,
):

    date, timestamp = obs_stamp

    if plots == "True":
        Path(f"{plt_dir}/{date}/{timestamp}").mkdir(parents=True, exist_ok=True)

    channel_map = {}

    ref_file = f"{ali_dir}/{date}/{timestamp}/rf0XX_S07XX_{timestamp}_aligned.npz"

    if Path(ref_file).is_file():
        ref_file = f"{ali_dir}/{date}/{timestamp}/rf0XX_S07XX_{timestamp}_aligned.npz"
    elif Path(
        ref_file=f"{ali_dir}/{date}/{timestamp}/rf0XX_S36XX_{timestamp}_aligned.npz"
    ).is_file():
        ref_file = f"{ali_dir}/{date}/{timestamp}/rf0XX_S36XX_{timestamp}_aligned.npz"
    else:
        ref_file = f"{ali_dir}/{date}/{timestamp}/rf0XX_S06XX_{timestamp}_aligned.npz"

    chrono_file = f"{chrono_dir}/{timestamp}.json"

    try:
        Path(ref_file).is_file()

        with open(chrono_file) as chrono:
            chrono_ephem = json.load(chrono)

            if chrono_ephem != []:

                norad_list = [
                    chrono_ephem[s]["sat_id"][0] for s in range(len(chrono_ephem))
                ]

                if norad_list != []:

                    for sat in norad_list:

                        sat_data = good_chans(
                            ref_file,
                            chrono_file,
                            sat,
                            sat_thresh,
                            noi_thresh,
                            pow_thresh,
                            occ_thresh,
                            date,
                            timestamp,
                            plots,
                        )

                        if sat_data != 0:
                            sat_chan = sat_data

                            channel_map[f"{sat}"] = sat_chan

    except Exception as e:
        print(e)
        # Exception message is forwarded from ../decode_rf_data/rf_data.py

    # Save channel map
    Path(f"{out_dir}/window_maps").mkdir(parents=True, exist_ok=True)
    with open(f"{out_dir}/window_maps/{timestamp}.json", "w") as f:
        json.dump(channel_map, f, indent=4)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
        Project a sat pass onto a healpix map using ephemeris data
        """
    )

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
        default="./../../outputs/sat_channels/",
        help="Output directory. Default=./../../outputs/sat_channels/",
    )
    parser.add_argument(
        "--plt_dir",
        metavar="\b",
        default="./../../outputs/sat_channels/window_plots",
        help="Output directory. Default=./../../outputs/sat_channels/window_plots/",
    )
    parser.add_argument(
        "--ali_dir",
        metavar="\b",
        default="./../../outputs/align_data/",
        help="Output directory. Default=./../../outputs/align_data/",
    )
    parser.add_argument(
        "--chrono_dir",
        metavar="\b",
        default="./../../outputs/sat_ephemeris/chrono_json",
        help="Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/",
    )
    parser.add_argument(
        "--noi_thresh",
        metavar="\b",
        type=int,
        default=3,
        help="Noise Threshold: Multiples of MAD. Default=3.",
    )
    parser.add_argument(
        "--sat_thresh",
        metavar="\b",
        type=int,
        default=1,
        help="1 σ threshold to detect sats Default=1.",
    )
    parser.add_argument(
        "--pow_thresh",
        metavar="\b",
        type=int,
        default=15,
        help="Power Threshold to detect sats. Default=15 dB.",
    )
    parser.add_argument(
        "--occ_thresh",
        metavar="\b",
        type=int,
        default=0.80,
        help="Occupation Threshold of sat in window. Default=0.80",
    )
    parser.add_argument(
        "--plots",
        metavar="\b",
        default=False,
        help="If True, create a gazzillion plots for each sat pass. Default = False",
    )

    args = parser.parse_args()

    chrono_dir = args.chrono_dir
    start_date = args.start_date
    stop_date = args.stop_date
    out_dir = args.out_dir
    plt_dir = args.plt_dir
    ali_dir = args.ali_dir
    noi_thresh = args.noi_thresh
    sat_thresh = args.sat_thresh
    pow_thresh = args.pow_thresh
    occ_thresh = args.occ_thresh
    plots = args.plots

    ref_names = ["rf0XX", "rf0YY", "rf1XX", "rf1YY"]

    # Save logs
    Path(f"{out_dir}/window_maps").mkdir(parents=True, exist_ok=True)

    # Help traverse all 30 min obs b/w start & stop
    dates, date_time = time_tree(start_date, stop_date)
    obs_list = [
        [dates[d], date_time[d][dt]]
        for d in range(len(dates))
        for dt in range(len(date_time[d]))
    ]

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(
            window_chan_map,
            repeat(ali_dir),
            repeat(chrono_dir),
            repeat(sat_thresh),
            repeat(noi_thresh),
            repeat(pow_thresh),
            repeat(occ_thresh),
            repeat(out_dir),
            repeat(plots),
            repeat(plt_dir),
            obs_list,
        )
