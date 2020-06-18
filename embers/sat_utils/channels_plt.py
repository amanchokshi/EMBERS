# Force matplotlib to not use X-Server backend
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
import seaborn as sns
import numpy as np


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
