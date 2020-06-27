"""
Satellite Channels
------------------

A set of tools to determine the transmission channel of various satellites in 30 minute
observation by using rf data in conjunctin with chronological satellite ephemeris data.

"""

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


def plt_waterfall_pass(
    power, sat_id, start, stop, date, cmap, chs=None, good_ch=None, out_dir=None
):
    """Plot waterfall with sat window and occupied channels
    
    Args:
        power:          RF power array
        sat_id:         Norad cat ID
        start:          Start of epehm for sat_id
        stop:           Stop of ephem of sat_id
        chs:            Occupied channels [list]
        good_ch:        Most probable channel
        date:           Date of observation
    """

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

    # horizontal highlight of ephem
    ax.axhspan(start, stop, alpha=0.1, color="white")

    if chs != None:
        # vertical highlight of channel
        for ch in chs:
            ax.axvspan(ch - 1.0, ch + 0.6, alpha=0.2, color="white")

    if good_ch != None:
        # highlight good channel
        ax.axvspan(
            good_ch - 1.0,
            good_ch + 0.6,
            alpha=0.6,
            edgecolor="#a7ff83",
            facecolor=None,
            fill=False,
        )

    # plt.show()
    if chs != None:
        plt.savefig(f"{out_dir}/{date}_{sat_id}_{chs}_waterfall.png")
    else:
        plt.savefig(f"{out_dir}/{date}_{sat_id}_waterfall.png")

    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)


def plt_channel_basic(
    out_dir,
    times,
    channel_power,
    chan_num,
    min_s,
    max_s,
    noise_threshold,
    arbitrary_threshold,
    sat_id,
    date,
):

    """Plot power in channel, with various thresholds
    
    Args:
        times:          Time array
        channel_power:  Power in channel
        chan_num:       Channel Number
        min_s:          Minimum signal in channel_power
        max_s:          Maximum signal in channel_power
        noise_threshold: Noise Threshold (n*MAD)
        arbitrary_threshold: Arbitrary threshold used to only select bright passes
        sat_id:         Norad Cat ID
        date:           Date of observation
        """

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
    plt.fill_between(times, channel_power, color="#db3751", alpha=0.7)

    plt.axhline(
        arbitrary_threshold,
        alpha=1.0,
        linestyle="-",
        linewidth=2,
        color="#fba95f",
        label=f"Arbitrary Cut: {arbitrary_threshold} dBm",
    )
    plt.axhspan(-1, arbitrary_threshold, color="#fba95f", alpha=0.4)

    plt.axhline(
        noise_threshold,
        linestyle="-",
        linewidth=2,
        color="#5cb7a9",
        label=f"Noise Cut: {noise_threshold:.2f} dBm",
    )
    plt.axhspan(-1, noise_threshold, color="#5cb7a9", alpha=0.4)

    plt.ylim([min_s - 1, max_s + 1])
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
    plt.savefig(f"{out_dir}/{date}_{sat_id}_{chan_num}_channel.png")
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)


def sat_plot(out_dir, ids, norad_id, alt, az, num_passes, date, name):
    """Plots satellite passes
    
    Args:
        alt: list of altitude values
        az: list of azimuth values
        num_passes: Number of satellite passes
    """

    # Set up the polar plot.
    plt.style.use("seaborn")
    plt.style.use("dark_background")
    figure = plt.figure(figsize=(7, 6))
    ax = figure.add_subplot(111, polar=True)
    ax.set_ylim(90, 0)
    ax.set_rgrids([0, 30, 60, 90], angle=22)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_title(f"{num_passes} satellite passes in {date} [{norad_id}] window", y=1.05)
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
        title="Norad SatID",
    )
    leg.get_frame().set_facecolor("grey")
    leg.get_frame().set_alpha(0.4)
    for l in leg.legendHandles:
        l.set_alpha(1)

    plt.tight_layout()
    plt.savefig(f"{out_dir}/{date}_{norad_id}_{name}.png")
    plt.close()
    plt.rcParams.update(plt.rcParamsDefault)


def good_chans(
    ref_file,
    chrono_file,
    sat_id,
    sat_thresh,
    noi_thresh,
    pow_thresh,
    occ_thresh,
    date,
    timestamp,
    plots,
):

    """Create power, alt, az arrays at constant cadence"""

    spec, _ = spectral()

    power, times = read_aligned(ali_file=ref_file)

    # Scale noise floor to zero and determine noise threshold
    power, noise_threshold = noise_floor(sat_thresh, noi_thresh, power)

    with open(chrono_file) as chrono:
        chrono_ephem = json.load(chrono)

        norad_list = [chrono_ephem[s]["sat_id"][0] for s in range(len(chrono_ephem))]

        norad_index = norad_list.index(sat_id)

        norad_ephem = chrono_ephem[norad_index]

        rise_ephem = norad_ephem["time_array"][0]
        set_ephem = norad_ephem["time_array"][-1]

        intvl = time_filter(rise_ephem, set_ephem, np.asarray(times))

        if intvl != None:

            w_start, w_stop = intvl

            window_len = w_stop - w_start + 1

            # Slice [crop] the power/times arrays to the times of sat pass
            power_c = power[w_start : w_stop + 1, :]
            times_c = times[w_start : w_stop + 1]

            possible_chans = []
            occu_list = []

            # Loop over every channel
            for s_chan in range(len(power_c[0])):

                channel_power = power_c[:, s_chan]

                # Percentage of signal occupancy above noise threshold
                window_occupancy = (np.where(channel_power >= noise_threshold))[
                    0
                ].size / window_len
                min_s = np.amin(channel_power)

                # Power threshold below which satellites aren't counted
                # Only continue if there is signal for more than 80% of satellite pass
                if (
                    max(channel_power) >= pow_thresh
                    and occ_thresh <= window_occupancy < 1.00
                ):

                    if times[0] == times_c[0]:
                        if (
                            all(p < noise_threshold for p in channel_power[-11:-1])
                        ) is True:
                            occu_list.append(window_occupancy)
                            possible_chans.append(s_chan)
                            # print(f'{sat_id}: {s_chan} {window_occupancy}')

                            if plots == "True":
                                plt_channel_basic(
                                    f"{plt_dir}/{date}/{timestamp}",
                                    times_c,
                                    channel_power,
                                    s_chan,
                                    min_s,
                                    32,
                                    noise_threshold,
                                    pow_thresh,
                                    sat_id,
                                    timestamp,
                                )

                    elif times[-1] == times_c[-1]:
                        if (
                            all(p < noise_threshold for p in channel_power[:10])
                        ) is True:
                            occu_list.append(window_occupancy)
                            possible_chans.append(s_chan)
                            # print(f'{sat_id}: {s_chan} {window_occupancy}')

                            if plots == "True":
                                plt_channel_basic(
                                    f"{plt_dir}/{date}/{timestamp}",
                                    times_c,
                                    channel_power,
                                    s_chan,
                                    min_s,
                                    32,
                                    noise_threshold,
                                    pow_thresh,
                                    sat_id,
                                    timestamp,
                                )

                    else:
                        if (
                            all(p < noise_threshold for p in channel_power[:10])
                            and all(p < noise_threshold for p in channel_power[-11:-1])
                        ) is True:
                            occu_list.append(window_occupancy)
                            possible_chans.append(s_chan)
                            # print(f'{sat_id}: {s_chan} {window_occupancy}')

                            if plots == "True":
                                plt_channel_basic(
                                    f"{plt_dir}/{date}/{timestamp}",
                                    times_c,
                                    channel_power,
                                    s_chan,
                                    min_s,
                                    32,
                                    noise_threshold,
                                    pow_thresh,
                                    sat_id,
                                    timestamp,
                                )

            # If channels are identified in the 30 min obs
            n_chans = len(possible_chans)

            if n_chans > 0:

                # The most lightly channel is one with the highest occupation
                good_chan = possible_chans[occu_list.index(max(occu_list))]

                if plots == "True":
                    plt_waterfall_pass(
                        power,
                        sat_id,
                        w_start,
                        w_stop,
                        timestamp,
                        spec,
                        chs=possible_chans,
                        good_ch=good_chan,
                        out_dir=f"{plt_dir}/{date}/{timestamp}",
                    )

                return good_chan

            else:
                if plots == "True":
                    plt_waterfall_pass(
                        power,
                        sat_id,
                        w_start,
                        w_stop,
                        timestamp,
                        spec,
                        chs=None,
                        good_ch=None,
                        out_dir=f"{plt_dir}/{date}/{timestamp}",
                    )
                return 0

        else:
            return 0


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
