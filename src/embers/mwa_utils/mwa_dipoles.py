"""
MWA Dipoles
-----------
Tools to download MWA metafits files and finding dead dipoles

"""

import json
import time
from pathlib import Path

import matplotlib as mpl
import numpy as np
import wget
from astropy.io import fits
from matplotlib import pylab as pl
from matplotlib import pyplot as plt

mpl.use("Agg")


def download_metafits(num_files, wait, out_dir):
    """
    Download metafits files from `mwatelescope.org <http://mwatelescope.org/>`_

    This function requires a list of obsids which it reads from :samp:`out_dir`.
    Before running functions in this module, run :mod:`~embers.mwa_utils.mwa_pointings` to
    created the required obsid files.

    :param num_files: The number of metafits files to download. Usually, 20 is sufficent :class:`~int`
    :param wait: Time to sleep between downloads so as not to overload servers. Default=29
    :param out_dir: Path to output directory where the metafits files will be saved :class:`~str`

    :returns:
        Metafits files are saved to the :samp:`out_dir`

    """

    print("Downloading MWA metafits files")
    print("Due to download limits, this will take a while")
    t = wait * num_files
    m, _ = divmod(t, 60)
    h, m = divmod(m, 60)
    print(f"ETA: Approximately {h:d}H:{m:02d}M")
    metafits_dir = Path(f"{out_dir}/mwa_metafits")
    metafits_dir.mkdir(parents=True, exist_ok=True)

    cerberus_url = "http://ws.mwatelescope.org/metadata/fits?obs_id="
    with open(f"{out_dir}/mwa_pointings.json") as gps:

        # list of all gpstimes [obsids]
        gps_times = np.array(json.load(gps)["start_gps"])

        # select 20 equally spaced ids
        idx = np.round(np.linspace(0, len(gps_times) - 1, num_files)).astype(int)
        obs_ids = gps_times[idx]

        for obs_id in obs_ids:

            time.sleep(wait)
            print(f"\nDownloading {obs_id} metafits")
            wget.download(
                f"{cerberus_url}{obs_id}", f"{metafits_dir}/{obs_id}.metafits"
            )

        print("\nMetafits download complete")


def find_flags(out_dir):
    """
    Read metafits files and determine which dipoles are flagged

    :param out_dir: Path to root of output directory where the metafits files are saved :class:`~str`

    :returns:
        A plot of flagged dipoles in each tile

    """

    # list of all tiles used in this experiment
    # Naming convetion used in metafit files
    tiles = [
        "HexS6",
        "HexS7",
        "HexS8",
        "HexS9",
        "HexS10",
        "HexS12",
        "HexS29",
        "HexS30",
        "HexS31",
        "HexS32",
        "HexS33",
        "HexS34",
        "HexS35",
        "HexS36",
    ]

    # dictionary to save flag info.
    # 0 indicates no dead dipole
    # 1-16 dipole index
    flags = {}
    flags["obsid"] = []
    for tile in tiles:
        for pol in ["X", "Y"]:
            flags[f"{tile}{pol}"] = []

    meta_files = [item for item in Path(f"{out_dir}/mwa_metafits/").glob("*.metafits")]
    for m in meta_files:

        hdu = fits.open(m)
        mode = hdu[0].header["MODE"]
        obsid = hdu[0].header["GPSTIME"]
        delays = hdu[1].data["Delays"]

        tile_names = hdu[1].data["TileName"]
        delays = hdu[1].data["Delays"]
        pols = hdu[1].data["Pol"]

        if mode == "HW_LFILES":

            flags["obsid"].append(int(obsid))

            for t in tiles:
                idx = np.where(tile_names == t)
                for i in range(2):
                    idx_p = idx[0][i]

                    t_name = tile_names[idx_p]
                    t_pol = pols[idx_p]
                    t_del = delays[idx_p]
                    t_flag = np.where(t_del == 32)[0] + 1

                    if t_flag.size != 0:
                        flags[f"{t_name}{t_pol}"].append(int(t_flag[0]))
                    else:
                        flags[f"{t_name}{t_pol}"].append(0)

        hdu.close()

    keys = list(flags.keys())

    n = len(keys) - 1
    colors = pl.cm.Spectral(np.linspace(0, 1, n))

    plt.style.use("seaborn")

    _, axs = plt.subplots(4, 7, figsize=(18, 9), sharex=True, sharey=True,)
    axs = axs.ravel()

    for i in range(n):
        axs[i].scatter(
            flags["obsid"],
            flags[keys[i + 1]],
            color=colors[i],
            linewidths=0.1,
            s=49,
            alpha=0.88,
            edgecolors="black",
            linewidth=0.6,
            label=keys[i + 1],
        )
        axs[i].set_ylim(-1, 17)
        axs[i].set_yticks([0, 4, 8, 12, 16])

        # these are matplotlib.patch.Patch properties
        props = dict(boxstyle="round", facecolor="wheat", alpha=0.2)

        # place a text box in upper left in axes coords
        axs[i].text(
            0.67,
            0.93,
            f"{keys[i+1]}",
            transform=axs[i].transAxes,
            fontsize=10,
            verticalalignment="top",
            bbox=props,
        )

    plt.tight_layout()
    plt.savefig(f"{out_dir}/flagged_dipoles.png")
    print(f"Dipole flagging plot saved to {out_dir}")


def mwa_flagged_dipoles(num_files, out_dir, wait=29):
    """Download metafits and find flagged dipoles

    :param num_files: The number of metafits files to download. Usually, 20 is sufficent :class:`~int`
    :param out_dir: Path to output directory where the metafits files will be saved :class:`~str`
    :param wait: Time to sleep between downloads so as not to overload servers. Default=29

    :returns:
        Metafits files and dipole flagging plot saved to :samp:`out_dir`

    """

    download_metafits(num_files, wait, out_dir)
    find_flags(out_dir)
