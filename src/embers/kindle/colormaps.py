"""
Colormaps
=========

Visualise custom colormaps used by embers.
Creates sample plot of ember colormaps
saved to :samp:`./embers_out/rf_tools/colormaps.png`

"""

import argparse
from pathlib import Path

from embers.rf_tools.colormaps import jade, plt_colormaps, spectral

_spec, _spec_r = spectral()
_jade, _jade_r = jade()


_parser = argparse.ArgumentParser(
    description="""
    Visualization of custom ember colormaps
    """
)

_parser.add_argument(
    "--out_dir",
    metavar="\b",
    default="./embers_out/rf_tools",
    help="Dir where colormap sample plot is saved. Default=./embers_out/rf_tools",
)

_args = _parser.parse_args()
_out_dir = Path(_args.out_dir)

# Make outdir if it doesn't exist
_out_dir.mkdir(parents=True, exist_ok=True)


def main():
    """Execute colormaps from terminal."""

    print(f"Plot of embers colormaps saved to ./{_out_dir}/colormaps.png")
    plt_colormaps(_spec, _spec_r, _jade, _jade_r, _out_dir)
