import sys
import json
import argparse
import matplotlib
import numpy as np
import healpy as hp

matplotlib.use("Agg")
from pathlib import Path
import concurrent.futures
import scipy.optimize as opt
from null_test import rotate
import matplotlib.pyplot as plt
from scipy.stats import chisquare
from scipy.stats import median_absolute_deviation as mad

sys.path.append("../sat_ephemeris")
from sat_ids import norad_ids

sys.path.append("../sat_channels")
from sat_channels import time_tree, time_filter

# Custom spectral colormap
sys.path.append("../decode_rf_data")
from rf_data import tile_names
from colormap import spectral

cmap = spectral()


def check_pointing(timestamp, point_0, point_2, point_4, point_41):
    """Check if timestamp is at pointing 0, 2, 4, 41"""
    if timestamp in point_0:
        point = 0
    elif timestamp in point_2:
        point = 2
    elif timestamp in point_4:
        point = 4
    else:
        point = 41

    return point


def read_aligned(ali_file=None):
    """Read aligned data npz file"""
    paired_data = np.load(ali_file, allow_pickle=True)

    ref_p = paired_data["ref_p_aligned"]
    tile_p = paired_data["tile_p_aligned"]
    times = paired_data["time_array"]

    return [ref_p, tile_p, times]


def noise_floor(sat_thresh, noi_thresh, data=None):
    """Computes the noise floor of the data """

    # compute the standard deviation of data, and use it to identify occupied channels
    σ = np.std(data)

    # Any channel with a max power >= σ has a satellite
    sat_cut = sat_thresh * σ
    chans_pow_max = np.amax(data, axis=0)

    # Exclude the channels with sats, to only have noise data
    noise_chans = np.where(chans_pow_max < sat_cut)[0]
    noise_data = data[:, noise_chans]

    # noise median, noise mad, noise threshold = μ + 3*σ
    μ_noise = np.median(noise_data)
    σ_noise = mad(noise_data, axis=None)
    noise_threshold = μ_noise + noi_thresh * σ_noise

    return (data, noise_threshold)


# chisquared minimization to best fit map to data
def fit_gain(map_data=None, fee=None):
    """Fit the beam model to the measured data using
    chisquared minimization"""

    bad_values = np.isnan(map_data)
    map_data = map_data[~bad_values]
    fee = fee[~bad_values]

    def chisqfunc(gain):
        model = fee + gain
        chisq = sum((map_data - model) ** 2)
        return chisq

    x0 = np.array([0])

    result = opt.minimize(chisqfunc, x0)

    return result.x


def fit_test(map_data=None, fee=None):
    """chi-squared test for fit between
    model and data"""

    bad_values = np.isnan(map_data)
    map_data = map_data[~bad_values]
    fee = fee[~bad_values]

    # map_data = np.asarray(map_data) + 120
    # fee = np.asarray(fee) + 120
    map_data = np.asarray(map_data) - np.nanmin(fee) + 20
    fee = np.asarray(fee) - np.nanmin(fee) + 20

    _, pvalue = chisquare(map_data, f_exp=fee)

    return pvalue


def plt_channel(
    out_dir,
    times,
    ref_c,
    tile_c,
    ref_noise,
    tile_noise,
    chan_num,
    sat_id,
    pointing,
    date,
    point,
):

    """Plot power in channel, with various thresholds
    
    Args:
        times:          Time array
        channel_power:  Power in channel
        chan_num:       Channel Number
        min_s:          Minimum signal in channel_power
        max_s:          Maximum signal in channel_power
        noise_threshold: Noise Threshold (n*MAD)
        sat_id:         Norad Cat ID
        date:           Date of observation
        """

    plt.style.use("seaborn")

    ylims = [min(ref_c), max(ref_c), min(tile_c), max(tile_c)]
    min_s = min(ylims)
    max_s = max(ylims)

    fig = plt.figure(figsize=(8, 6))
    fig.suptitle(f"Satellite [{sat_id}] Pass @ {date} in Channel: [{chan_num}]", y=0.98)

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(
        times,
        ref_c,
        linestyle="-",
        linewidth=2,
        alpha=1.0,
        color="#729d39",
        label="ref",
    )
    ax1.fill_between(times, y1=ref_c, y2=-120, color="#729d39", alpha=0.7)
    ax1.axhline(
        ref_noise,
        linestyle="-",
        linewidth=2,
        color="#36622b",
        label=f"Ref Cut: {ref_noise:.2f} dBm",
    )
    ax1.axhspan(-120, ref_noise, color="#36622b", alpha=0.7)
    ax1.set_ylim([min(ref_c) - 1, max(ref_c) + 1])
    ax1.set_ylabel("Power [dBm]")
    ax1.set_xticklabels([])

    leg = ax1.legend(frameon=True)
    leg.get_frame().set_facecolor("grey")
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(
        times,
        tile_c,
        linestyle="-",
        linewidth=2,
        alpha=1.0,
        color="#ff5656",
        label="tile",
    )
    ax2.fill_between(times, y1=tile_c, y2=-120, color="#ff5656", alpha=0.7)
    ax2.axhline(
        tile_noise,
        linestyle="-",
        linewidth=2,
        color="#970747",
        label=f"Tile Cut: {tile_noise:.2f} dBm",
    )
    ax2.axhspan(-120, tile_noise, color="#970747", alpha=0.7)
    ax2.set_ylim([min(tile_c) - 1, max(tile_c) + 1])
    ax2.set_ylabel("Power [dBm]")
    ax2.set_xlabel("Time [s]")

    leg = ax2.legend(frameon=True)
    leg.get_frame().set_facecolor("grey")
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)

    plt.tight_layout()
    plt.subplots_adjust(top=0.94)
    plt.savefig(f"{out_dir}/{date}_{sat_id}_{chan_num}_channel.png")
    plt.close()


def plt_fee_fit(
    t, mwa_fee_pass, mwa_pass_fit, mwa_pass_fit_raw, out_dir, point, timestamp, sat
):

    pval = fit_test(map_data=mwa_pass_fit, fee=mwa_fee_pass)

    plt.style.use("seaborn")

    fig = plt.figure(figsize=(8, 6))
    ax1 = fig.add_subplot(1, 1, 1)

    ax1.scatter(
        t, mwa_fee_pass, color="#c70039", alpha=0.6, marker=".", label="fee_slice"
    )

    where_nan = np.isnan(mwa_pass_fit)

    ax1.scatter(
        t, mwa_pass_fit_raw, color="#7da87b", alpha=0.9, marker=".", label="RFE raw"
    )
    ax1.set_ylim(-50, 2)

    ax1.scatter(
        t, mwa_pass_fit, color="#4f8a8b", alpha=0.9, marker=".", label="RFE cali"
    )
    ax1.set_ylim(-50, 2)

    leg = ax1.legend(frameon=True)
    leg.get_frame().set_facecolor("grey")
    leg.get_frame().set_alpha(0.2)
    for l in leg.legendHandles:
        l.set_alpha(1)

    plt.title(f"Goodness of Fit: {pval} ")
    plt.tight_layout()
    plt.savefig(f"{out_dir}/{point}/{timestamp}_{sat}.png")
    plt.close()


def power_ephem(
    ref, tile, ali_file, chrono_file, sat_id, sat_chan, point, pow_thresh, timestamp
):

    """Create power, alt, az arrays at constant cadence"""

    # Read .npz aligned file
    ref_p, tile_p, times = read_aligned(ali_file=ali_file)

    # Scale noise floor to zero and determine noise threshold
    ref_p, ref_noise = noise_floor(sat_thresh, noi_thresh, ref_p)
    tile_p, tile_noise = noise_floor(sat_thresh, noi_thresh, tile_p)

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

            # Slice [crop] the ref/tile/times arrays to the times of sat pass and extract sat_chan
            ref_c = ref_p[w_start : w_stop + 1, sat_chan]
            tile_c = tile_p[w_start : w_stop + 1, sat_chan]
            times_c = times[w_start : w_stop + 1]

            alt = np.asarray(norad_ephem["sat_alt"])
            az = np.asarray(norad_ephem["sat_az"])

            if (np.nanmax(ref_c - ref_noise) >= pow_thresh) and (
                np.nanmax(tile_c - tile_noise) >= pow_thresh
            ):

                # Apply noise criteria. In the window, where are ref_power and tile power
                # above their respective thresholds?
                if np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0].size != 0:
                    good_ref = ref_c[
                        np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]
                    ]
                    good_tile = tile_c[
                        np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]
                    ]
                    good_alt = alt[
                        np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]
                    ]
                    good_az = az[
                        np.where((ref_c >= ref_noise) & (tile_c >= tile_noise))[0]
                    ]

                    if plots == "True":
                        plt_channel(
                            f"{out_dir}/pass_plots/{tile}_{ref}/{point}",
                            times_c,
                            ref_c,
                            tile_c,
                            ref_noise,
                            tile_noise,
                            sat_chan,
                            sat_id,
                            point,
                            timestamp,
                            point,
                        )

                    return [good_ref, good_tile, good_alt, good_az, times_c]

                else:
                    return 0

            else:
                return 0

        else:
            return 0


def project_tile_healpix(tile_pair):

    ref, tile = tile_pair

    pointings = ["0", "2", "4", "41"]

    if plots == "True":
        Path(f"{out_dir}/pass_plots/{tile}_{ref}/0").mkdir(parents=True, exist_ok=True)
        Path(f"{out_dir}/pass_plots/{tile}_{ref}/2").mkdir(parents=True, exist_ok=True)
        Path(f"{out_dir}/pass_plots/{tile}_{ref}/4").mkdir(parents=True, exist_ok=True)
        Path(f"{out_dir}/pass_plots/{tile}_{ref}/41").mkdir(parents=True, exist_ok=True)

        Path(f"{out_dir}/fit_plots/{tile}_{ref}/0").mkdir(parents=True, exist_ok=True)
        Path(f"{out_dir}/fit_plots/{tile}_{ref}/2").mkdir(parents=True, exist_ok=True)
        Path(f"{out_dir}/fit_plots/{tile}_{ref}/4").mkdir(parents=True, exist_ok=True)
        Path(f"{out_dir}/fit_plots/{tile}_{ref}/41").mkdir(parents=True, exist_ok=True)

    # Initialize an empty dictionary for tile data
    # The map is list of length 12288 of empty lists to append pixel values to
    # keep track of which satellites contributed which data
    tile_data = {
        "mwa_maps": {
            p: [[] for pixel in range(hp.nside2npix(nside))] for p in pointings
        },
        "ref_maps": {
            p: [[] for pixel in range(hp.nside2npix(nside))] for p in pointings
        },
        "tile_maps": {
            p: [[] for pixel in range(hp.nside2npix(nside))] for p in pointings
        },
        "sat_map": {
            p: [[] for pixel in range(hp.nside2npix(nside))] for p in pointings
        },
        "times": {p: [[] for pixel in range(hp.nside2npix(nside))] for p in pointings},
    }

    # Load reference FEE model
    # Rotate the fee models by -pi/2 to move model from spherical (E=0) to Alt/Az (N=0)
    ref_fee_model = np.load(ref_model, allow_pickle=True)

    #    # Tile S33YY has it's 9th dipole flagged
    #    if tile is 'S33YY':
    #        fee_m = np.load(fee_map_flagged, allow_pickle=True)
    #    else:
    #        fee_m = np.load(fee_map, allow_pickle=True)
    fee_m = np.load(fee_map, allow_pickle=True)

    if "XX" in tile:
        ref_fee = ref_fee_model["XX"]
        rotated_fee = rotate(nside, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee)
    else:
        ref_fee = ref_fee_model["YY"]
        rotated_fee = rotate(nside, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee)

    for day in range(len(dates)):

        for window in range(len(date_time[day])):
            timestamp = date_time[day][window]

            # Check if at timestamp, reciever was pointed to 0,2,4 gridpointing
            if (
                (timestamp in point_0)
                or (timestamp in point_2)
                or (timestamp in point_4)
                or (timestamp in point_41)
            ):

                # pointing at timestamp
                point = check_pointing(timestamp, point_0, point_2, point_4, point_41)

                if "XX" in tile:
                    mwa_fee = fee_m[str(point)][0]
                else:
                    mwa_fee = fee_m[str(point)][1]

                ali_file = Path(
                    f"{align_dir}/{dates[day]}/{timestamp}/{ref}_{tile}_{timestamp}_aligned.npz"
                )

                # check if file exists
                if ali_file.is_file():

                    # Chrono and map Ephemeris file
                    chrono_file = Path(f"{chrono_dir}/{timestamp}.json")
                    channel_map = f"{map_dir}/{timestamp}.json"

                    with open(chrono_file) as chrono:
                        chrono_ephem = json.load(chrono)

                        if chrono_ephem != []:

                            norad_list = [
                                chrono_ephem[s]["sat_id"][0]
                                for s in range(len(chrono_ephem))
                            ]

                            if norad_list != []:

                                with open(channel_map) as ch_map:
                                    chan_map = json.load(ch_map)

                                    chan_sat_ids = [
                                        int(i) for i in list(chan_map.keys())
                                    ]

                                    for sat in chan_sat_ids:

                                        chan = chan_map[f"{sat}"]

                                        sat_data = power_ephem(
                                            ref,
                                            tile,
                                            ali_file,
                                            chrono_file,
                                            sat,
                                            chan,
                                            point,
                                            pow_thresh,
                                            timestamp,
                                        )

                                        if sat_data != 0:

                                            (
                                                ref_power,
                                                tile_power,
                                                alt,
                                                az,
                                                times,
                                            ) = sat_data

                                            # Altitude is in deg while az is in radians
                                            # convert alt to radians
                                            # za - zenith angle
                                            alt = np.radians(alt)
                                            za = np.pi / 2 - alt
                                            az = np.asarray(az)

                                            # Now convert to healpix coordinates
                                            # healpix_index = hp.ang2pix(nside,θ, ɸ)
                                            healpix_index = hp.ang2pix(nside, za, az)

                                            # multiple data points fall within a single healpix pixel
                                            # find the unique pixels
                                            u = np.unique(healpix_index)

                                            ref_pass = np.array(
                                                [
                                                    np.nanmean(
                                                        ref_power[
                                                            np.where(
                                                                healpix_index == i
                                                            )[0]
                                                        ]
                                                    )
                                                    for i in u
                                                ]
                                            )
                                            tile_pass = np.array(
                                                [
                                                    np.nanmean(
                                                        tile_power[
                                                            np.where(
                                                                healpix_index == i
                                                            )[0]
                                                        ]
                                                    )
                                                    for i in u
                                                ]
                                            )
                                            times_pass = np.array(
                                                [
                                                    np.mean(
                                                        times[
                                                            np.where(healpix_index == i)
                                                        ][0]
                                                    )
                                                    for i in u
                                                ]
                                            )
                                            ref_fee_pass = np.array(
                                                [rotated_fee[i] for i in u]
                                            )
                                            mwa_fee_pass = np.array(
                                                [mwa_fee[i] for i in u]
                                            )

                                            # implement RFE gain corrections here
                                            rfe_polyfit = np.load(rfe_gain)
                                            gain_cal = np.poly1d(rfe_polyfit)
                                            # rfe_thresh = gain_cal.roots[0]
                                            rfe_thresh = max(gain_cal.roots)
                                            tile_pass_rfe = [
                                                i + gain_cal(i)
                                                if i >= rfe_thresh
                                                else i
                                                for i in tile_pass
                                            ]

                                            # magic here
                                            # the beam shape finally emerges
                                            # mwa_pass = np.array(tile_pass) - np.array(ref_pass) + np.array(ref_fee_pass)
                                            mwa_pass = (
                                                np.array(tile_pass_rfe)
                                                - np.array(ref_pass)
                                                + np.array(ref_fee_pass)
                                            )

                                            # fit the power level of the pass to the mwa_fee model using a single gain value
                                            offset = fit_gain(
                                                map_data=mwa_pass, fee=mwa_fee_pass
                                            )

                                            mwa_pass_fit = mwa_pass - offset[0]

                                            if mwa_pass_fit.size != 0:

                                                # determine how well the data fits the model with chi-square
                                                pval = fit_test(
                                                    map_data=mwa_pass_fit,
                                                    fee=mwa_fee_pass,
                                                )

                                                if plots == "True":
                                                    mwa_pass_raw = (
                                                        np.array(tile_pass)
                                                        - np.array(ref_pass)
                                                        + np.array(ref_fee_pass)
                                                    )
                                                    offset = fit_gain(
                                                        map_data=mwa_pass_raw,
                                                        fee=mwa_fee_pass,
                                                    )
                                                    mwa_pass_fit_raw = (
                                                        mwa_pass_raw - offset[0]
                                                    )

                                                    plt_fee_fit(
                                                        times_pass,
                                                        mwa_fee_pass,
                                                        mwa_pass_fit,
                                                        mwa_pass_fit_raw,
                                                        f"{out_dir}/fit_plots/{tile}_{ref}/",
                                                        point,
                                                        timestamp,
                                                        sat,
                                                    )

                                                # a goodness of fit threshold
                                                if pval >= 0.8:

                                                    # loop though all healpix pixels for the pass
                                                    for i in range(len(u)):

                                                        tile_data["mwa_maps"][
                                                            f"{point}"
                                                        ][u[i]].append(mwa_pass_fit[i])
                                                        tile_data["ref_maps"][
                                                            f"{point}"
                                                        ][u[i]].append(ref_pass[i])
                                                        tile_data["tile_maps"][
                                                            f"{point}"
                                                        ][u[i]].append(tile_pass[i])
                                                        tile_data["times"][f"{point}"][
                                                            u[i]
                                                        ].append(times_pass[i])
                                                        tile_data["sat_map"][
                                                            f"{point}"
                                                        ][u[i]].append(sat)

                else:
                    print(f"Missing {ref}_{tile}_{timestamp}_aligned.npz")
                    continue

    # Sort data by satellites

    # list of all possible satellites
    sat_ids = list(norad_ids.values())

    # create a dictionary, with the keys being sat_ids and the values being healpix maps of data from those sats
    # Fist level keys are pointings with dictionaty values
    # Second level keys are satellite ids with healpix map lists
    mwa_sat_data = {}
    ref_sat_data = {}
    tile_sat_data = {}
    time_sat_data = {}

    # loop over pointings for all maps
    for p in pointings:

        mwa_map = np.asarray(tile_data["mwa_maps"][p])
        tile_map = np.asarray(tile_data["tile_maps"][p])
        ref_map = np.asarray(tile_data["ref_maps"][p])
        sat_map = np.asarray(tile_data["sat_map"][p])
        time_map = np.asarray(tile_data["times"][p])

        mwa_sat_data[p] = {}
        ref_sat_data[p] = {}
        tile_sat_data[p] = {}
        time_sat_data[p] = {}

        # loop over all sats
        for s in sat_ids:

            mwa_sat_data[p][s] = []
            ref_sat_data[p][s] = []
            tile_sat_data[p][s] = []
            time_sat_data[p][s] = []

            # loop over every healpix pixel
            for i in range(len(ref_map)):

                # Find subset of data for each sat
                sat_idx = np.where(np.asarray(sat_map[i]) == s)
                mwa_sat_data[p][s].append((np.asarray(mwa_map[i])[sat_idx]).tolist())
                ref_sat_data[p][s].append((np.asarray(ref_map[i])[sat_idx]).tolist())
                tile_sat_data[p][s].append((np.asarray(tile_map[i])[sat_idx]).tolist())
                time_sat_data[p][s].append((np.asarray(time_map[i])[sat_idx]).tolist())

    # Save map arrays to npz file
    tile_sat_data = {
        "mwa_map": mwa_sat_data,
        "ref_map": ref_sat_data,
        "tile_map": tile_sat_data,
        "time_map": time_sat_data,
    }
    np.savez_compressed(f"{out_dir}/{tile}_{ref}_sat_maps.npz", **tile_sat_data)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
        Project a sat pass onto a healpix map using ephemeris data
        """
    )

    parser.add_argument(
        "--align_dir",
        metavar="\b",
        default="../../outputs/align_data/",
        help="Dir where agligned data date is saved. Default:../../outputs/align_data",
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./../../outputs/tile_maps/tile_maps_raw/",
        help="Output directory. Default=./../../outputs/tile_maps/tile_maps_raw/",
    )

    parser.add_argument(
        "--obs_point",
        metavar="\b",
        default="../../outputs/beam_pointings/obs_pointings.json",
        help="Observation pointing lists. Default=../../outputs/beam_pointings/obs_pointings.json",
    )

    parser.add_argument(
        "--chrono_dir",
        metavar="\b",
        default="./../../outputs/sat_ephemeris/chrono_json",
        help="Output directory. Default=./../../outputs/sat_ephemeris/chrono_json/",
    )

    parser.add_argument(
        "--map_dir",
        metavar="\b",
        default="../../outputs/sat_channels/window_maps/",
        help="Satellite channel map. Default=../../outputs/sat_channels/window_maps/",
    )

    parser.add_argument(
        "--ref_model",
        metavar="\b",
        default="../../outputs/reproject_ref/ref_dipole_models.npz",
        help="Healpix reference FEE model file. default=../../outputs/reproject_ref/ref_dipole_models.npz",
    )

    parser.add_argument(
        "--rfe_gain",
        metavar="\b",
        default="../../outputs/tile_maps/rfe_gain/rfe_gain_fit.npy",
        help="RF Explorer gain fit. default=../../outputs/tile_maps/rfe_gain/rfe_gain_fit.npy",
    )

    parser.add_argument(
        "--fee_map",
        metavar="\b",
        default="../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz",
        help="Healpix FEE map of mwa tile. default=../../outputs/tile_maps/FEE_maps/mwa_fee_beam.npz",
    )

    parser.add_argument(
        "--fee_map_flagged",
        metavar="\b",
        default="../../outputs/tile_maps/FEE_maps/mwa_fee_beam_9_flagged.npz",
        help="Healpix FEE map of mwa tile. default=../../outputs/tile_maps/FEE_maps/mwa_fee_beam_9_flagged.npz",
    )

    parser.add_argument(
        "--fit_thresh",
        metavar="\b",
        default=0.9,
        help="Goodness of fit threshold. 0.9 seems to only reject obvious outliers",
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
        help="σ threshold to detect sats Default=1.",
    )
    parser.add_argument(
        "--pow_thresh",
        metavar="\b",
        type=int,
        default=5,
        help="Power Threshold to detect sats. Default=10 dB.",
    )
    parser.add_argument(
        "--nside",
        metavar="\b",
        type=int,
        default=32,
        help="Healpix Nside. Default = 32",
    )
    parser.add_argument(
        "--plots",
        metavar="\b",
        default=False,
        help="If True, create a gazzillion plots for each sat pass. Default = False",
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

    args = parser.parse_args()

    start_date = args.start_date
    stop_date = args.stop_date
    obs_point = args.obs_point
    noi_thresh = args.noi_thresh
    sat_thresh = args.sat_thresh
    pow_thresh = args.pow_thresh
    fit_thresh = args.fit_thresh
    nside = args.nside
    plots = args.plots
    ref_model = args.ref_model
    rfe_gain = args.rfe_gain
    fee_map = args.fee_map
    fee_map_flagged = args.fee_map_flagged

    align_dir = Path(args.align_dir)
    chrono_dir = Path(args.chrono_dir)
    out_dir = Path(args.out_dir)
    map_dir = Path(args.map_dir)

    # Tile names
    refs = tile_names()[:4]
    tiles = tile_names()[4:]

    # All relevant tile pairs
    tile_pairs = []
    for ref in refs:
        if "XX" in ref:
            for tile in [t for t in tiles if "XX" in t]:
                tile_pairs.append([ref, tile])
        else:
            for tile in [t for t in tiles if "YY" in t]:
                tile_pairs.append([ref, tile])

    # Read observation pointing list
    with open(obs_point) as point:
        obs_p = json.load(point)
        point_0 = obs_p["point_0"]
        point_2 = obs_p["point_2"]
        point_4 = obs_p["point_4"]
        point_41 = obs_p["point_41"]

    # dates: list of days
    # date_time = list of 30 min observation windows
    dates, date_time = time_tree(start_date, stop_date)

    # Save logs
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    sys.stdout = open(f"{out_dir}/logs_{start_date}_{stop_date}.txt", "a")

    #    project_tile_healpix(tile_pairs[0])
    #    project_tile_healpix(['rf0YY', 'S33YY'])

    # Parallization magic happens here
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(project_tile_healpix, tile_pairs)
