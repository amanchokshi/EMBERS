"""Visualise custom colormaps used bt embers
.. code-block:: bash

   python -m embers.kindle.colormaps

:return: Saved plot showing included colormaps. By default save to ``./outputs/colormaps.png``
"""

import argparse
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from embers.decode_rfdata.colormaps import spectral, jade

spec, spec_r = spectral()
jade, jade_r = jade()


def waves_2d():
    """Creates 2d sine wave"""

    xx, yy = np.meshgrid(np.linspace(0, 3 * np.pi, 512), np.linspace(0, 3 * np.pi, 512))
    z = np.sin(xx) + np.sin(yy)
    return z


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""
        Visualization of custom ember colormaps
        """
    )

    parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./outputs",
        help="Dir where colormap sample plot is saved. Default=./",
    )

    args = parser.parse_args()
    out_dir = Path(args.out_dir)

    # Make outdir if it doesn't exist
    out_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(9, 8))

    ax1 = fig.add_subplot(221)
    ax1.set_title("Spectral colormap")
    im1 = ax1.imshow(waves_2d(), origin="lower", interpolation="none", cmap=spec)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im1, cax=cax, orientation="vertical")
    ax1.set_xticklabels([])

    ax2 = fig.add_subplot(222)
    ax2.set_title("Jade colormap")
    im2 = ax2.imshow(waves_2d(), origin="lower", interpolation="none", cmap=jade)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im2, cax=cax, orientation="vertical")
    ax2.set_xticklabels([])
    ax2.set_yticklabels([])

    ax3 = fig.add_subplot(223)
    ax3.set_title("Spectral_r colormap")
    im3 = ax3.imshow(waves_2d(), origin="lower", interpolation="none", cmap=spec_r)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im3, cax=cax, orientation="vertical")

    ax4 = fig.add_subplot(224)
    ax4.set_title("Jade_r colormap")
    im4 = ax4.imshow(waves_2d(), origin="lower", interpolation="none", cmap=jade_r)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(im4, cax=cax, orientation="vertical")
    ax4.set_yticklabels([])

    plt.tight_layout()
    plt.savefig(f"{out_dir}/colormaps.png")