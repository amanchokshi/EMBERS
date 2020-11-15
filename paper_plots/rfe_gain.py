"""
RFE GAIN.
---------

"""

import json
from pathlib import Path

import healpy as hp
import matplotlib
import numpy as np
from embers.sat_utils.sat_channels import time_tree
from embers.tile_maps.beam_utils import chisq_fit_gain, rotate_map
from embers.tile_maps.tile_maps import check_pointing, rf_apply_thresholds
from matplotlib import pyplot as plt
from embers.rf_tools.colormaps import spectral

_spec, _ = spectral()
colors = _spec([0.72, 0.30, 0.18])

matplotlib.use("Agg")

start_date = "2019-10-02"
stop_date = "2019-10-03"
tile_pair = ["rf0XX", "S06XX"]
sat_thresh = 1
noi_thresh = 3
pow_thresh = 5
ref_model = "../embers_out/tile_maps/ref_models/ref_dipole_models.npz"
fee_map = "../embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz"
rfe_cali = "../embers_out/tile_maps/rfe_calibration/rfe_gain_fit.npy"
nside = 32
obs_point_json = "../embers_out/mwa_utils/obs_pointings.json"
align_dir = "../embers_out/rf_tools/align_data"
chrono_dir = "../embers_out/sat_utils/ephem_chrono"
chan_map_dir = "../embers_out/sat_utils/sat_channels/window_maps"
out_dir = "../embers_out/paper_plots/rfe_gain"

Path(out_dir).mkdir(parents=True, exist_ok=True)


def plt_fee_fit(
    times, mwa_fee_pass, mwa_pass_fit_raw, mwa_pass_fit, out_dir, point, timestamp, sat
):
    """Plot data and model with goodness of fit p-value to visualize the degree of fit

    :param times: Time array
    :param mwa_fee_pass: MWA fee model slices according to satellite pass ephemeris
    :param mwa_pass_fit_raw: Satellite rf data from MWA tiles
    :param mwa_pass_fit: Satellite rf data from MWA tiles, fit to the fee model
    :param out_dir: Output directory where plot will be saved
    :param point: MWA sweet pointing of the observation
    :param timestamp:  Observation timestamp
    :param sat_id: Norad Cat ID

    :returns:
        - Plot comparing MWA fee model slice to satellite data, before and after gain power corrections

    """

    plt.style.use("seaborn")
    nice_fonts = {
        # Use LaTeX to write all text
        # "text.usetex": True,
        "font.family": "sans-serif",
        # Use 10pt font in plots, to match 10pt font in document
        "axes.labelsize": 10,
        "font.size": 10,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 6,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
    }

    plt.rcParams.update(nice_fonts)

    fig = plt.figure(figsize=(3.6, 2.4))
    ax1 = fig.add_subplot(1, 1, 1)

    times = (times - min(times))/60

    ax1.scatter(
        times, mwa_fee_pass, color=colors[0], alpha=0.6, marker="p", s=10, label=r"B$_{FEE}$"
        #  times, mwa_fee_pass, color="#c70039", alpha=0.6, marker="p", s=10, label=r"B$_{fee}$"
    )

    ax1.scatter(
        times,
        mwa_pass_fit_raw,
        color=colors[1],
        #  color="#7da87b",
        alpha=0.9,
        marker="p",
        s=10,
        label=r"B$_{def}$",
    )
    ax1.set_ylim(-50, 2)

    ax1.scatter(
        times,
        mwa_pass_fit,
        color=colors[2],
        #  color="#023e8a",
        alpha=0.9,
        marker="p",
        s=10,
        label=r"B$_{cali}$",
    )
    ax1.set_ylim(-50, 2)
    ax1.set_ylabel("Power [dBm]")
    ax1.set_xlabel("Time [minutes]")

    leg = ax1.legend(frameon=True, markerscale=2, loc="upper left")
    leg.get_frame().set_facecolor("white")
    for le in leg.legendHandles:
        le.set_alpha(1)

    plt.tight_layout()
    plt.savefig(f"{out_dir}/{timestamp}_{sat}.pdf", bbox_inches="tight")
    plt.close()


def project_tile_healpix(
    start_date,
    stop_date,
    tile_pair,
    sat_thresh,
    noi_thresh,
    pow_thresh,
    ref_model,
    fee_map,
    rfe_cali,
    nside,
    obs_point_json,
    align_dir,
    chrono_dir,
    chan_map_dir,
    out_dir,
    plots=True,
):
    """There be magic here. Project satellite RF data onto a sky healpix map.

    :param start_date: Start date in :samp:`YYYY-MM-DD-HH:MM` format
    :param stop_date: Stop date in :samp:`YYYY-MM-DD-HH:MM` format
    :param tile_pair: A pair of reference and MWA tile names. Ex: ["rf0XX", "S06XX"]
    :param sat_thresh: σ threshold to detect sats in the computation of rf data noise_floor. A good default is 1
    :param noi_thresh: Noise Threshold: Multiples of MAD. 3 is a good default
    :param pow_thresh: Peak power which must be exceeded for satellite pass to be considered
    :param ref_model: Path to reference feko model :samp:`.npz` file, output by :func:`~embers.tile_maps.ref_fee_healpix.ref_healpix_save`
    :param fee_map: Path to MWA fee model :samp:`.npz` file, output by :func:`~embers.mwa_utils.mwa_fee.mwa_fee_model`
    :param rfe_cali: Path to RFE gain calibration solution, output by :func:`~embers.tile_maps.tile_maps.rfe_collate_cali`
    :param nside: Healpix nside
    :param obs_point_json: Path to :samp:`obs_pointings.json` created by :func:`~embers.mwa_utils.mwa_pointings.obs_pointings`
    :param align_dir: Path to directory containing aligned rf data files, output from :func:`~embers.rf_tools.align_data.save_aligned`
    :param chrono_dir: Path to directory containing chronological ephemeris data output from :func:`~embers.sat_utils.chrono_ephem.save_chrono_ephem`
    :param chan_map_dir: Path to directory containing satellite frequency channel maps. Output from :func:`~embers.sat_utils.sat_channels.batch_window_map`
    :param out_dir: Output directory where rfe calibration data will be saved as a :samp:`json` file
    :param plots: If True, create a zillion diagnostic plots

    :returns:
        - Tile maps saved as :samp:`.npz` file to :samp:`out_dir`

    """

    ref, tile = tile_pair

    dates, timestamps = time_tree(start_date, stop_date)

    # Load reference FEE model
    # Rotate the fee models by -pi/2 to move model from spherical (E=0) to Alt/Az (N=0)
    ref_fee_model = np.load(ref_model, allow_pickle=True)

    fee_m = np.load(fee_map, allow_pickle=True)

    if "XX" in tile:
        ref_fee = ref_fee_model["XX"]
        rotated_fee = rotate_map(nside, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee)
    else:
        ref_fee = ref_fee_model["YY"]
        rotated_fee = rotate_map(nside, angle=-(1 * np.pi) / 2.0, healpix_array=ref_fee)

    for day in range(len(dates)):

        for window in range(len(timestamps[day])):
            timestamp = timestamps[day][window]

            # pointing at timestamp
            point = check_pointing(timestamp, obs_point_json)

            if point == 0:

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
                    channel_map = Path(f"{chan_map_dir}/{timestamp}.json")

                    with open(chrono_file) as chrono:
                        chrono_ephem = json.load(chrono)

                        if chrono_ephem != []:

                            norad_list = [
                                chrono_ephem[s]["sat_id"][0]
                                for s in range(len(chrono_ephem))
                            ]

                            if norad_list != []:

                                if channel_map.is_file():

                                    with open(channel_map) as ch_map:
                                        chan_map = json.load(ch_map)

                                        chan_sat_ids = [
                                            int(i) for i in list(chan_map.keys())
                                        ]

                                        for sat in chan_sat_ids:

                                            chan = chan_map[f"{sat}"]

                                            sat_data = rf_apply_thresholds(
                                                ali_file,
                                                chrono_file,
                                                sat,
                                                chan,
                                                sat_thresh,
                                                noi_thresh,
                                                pow_thresh,
                                                point,
                                                False,
                                                f"{out_dir}/tile_maps_raw",
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
                                                healpix_index = hp.ang2pix(
                                                    nside, za, az
                                                )

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
                                                                np.where(
                                                                    healpix_index == i
                                                                )
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
                                                # Read in RFE gain calibration solution
                                                rfe_polyfit = np.load(rfe_cali)
                                                gain_cal = np.poly1d(rfe_polyfit)

                                                # The max root of the polynomial is where we begin to apply the gain correction from
                                                rfe_thresh = max(gain_cal.roots)

                                                # When the power exceeds rfe_thresh, add to it using the rfe gain polynomial
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
                                                offset = chisq_fit_gain(
                                                    data=mwa_pass, model=mwa_fee_pass
                                                )

                                                mwa_pass_fit = mwa_pass - offset[0]

                                                if mwa_pass_fit.size != 0:

                                                    if plots is True:

                                                        if max(mwa_fee_pass) >= -10:

                                                            mwa_pass_raw = (
                                                                np.array(tile_pass)
                                                                - np.array(ref_pass)
                                                                + np.array(ref_fee_pass)
                                                            )
                                                            offset = chisq_fit_gain(
                                                                data=mwa_pass_raw,
                                                                model=mwa_fee_pass,
                                                            )
                                                            mwa_pass_fit_raw = (
                                                                mwa_pass_raw - offset[0]
                                                            )

                                                            plt_fee_fit(
                                                                times_pass,
                                                                mwa_fee_pass,
                                                                mwa_pass_fit_raw,
                                                                mwa_pass_fit,
                                                                out_dir,
                                                                point,
                                                                timestamp,
                                                                sat,
                                                            )


project_tile_healpix(
    start_date,
    stop_date,
    tile_pair,
    sat_thresh,
    noi_thresh,
    pow_thresh,
    ref_model,
    fee_map,
    rfe_cali,
    nside,
    obs_point_json,
    align_dir,
    chrono_dir,
    chan_map_dir,
    out_dir,
)
print(f"RFE GAIN plots saved to {out_dir}")
