"""
MWA FEE
-------
Tool to create XX & YY FEE (Fully Embedded Element) simulated beam maps.

"""

import argparse

from embers.mwa_utils.mwa_fee import mwa_fee_model


def main():
    """
    Create MWA Fully Embedded Element (FEE) beam models healpix maps at the given nside using the :func:`~embers.mwa_utils.mwa_fee.mwa_fee_model` function.

    .. code-block:: console

        $ mwa_fee --help
    """

    _parser = argparse.ArgumentParser(
        description="""
        Create MWA FEE beam models at multiple pointings with dipoles flagged.
        """
    )

    _parser.add_argument(
        "--nside", metavar="\b", default=32, type=int, help="Healpix nside. Default=32"
    )

    _parser.add_argument(
        "--pointings",
        metavar="\b",
        default="0, 2, 4, 41",
        type=str,
        help="List of MWA sweet spot pointings at which to evaluate the beam. Default='0, 2, 4, 41'",
    )

    _parser.add_argument(
        "--flags",
        metavar="\b",
        default="0",
        type=str,
        help="List of flagged dipoles with indices from 1 to 32. 1-16 are dipoles of XX pol while 17-32 are for YY. Ex: flags='1,17' represents the first dipole of the XX & YY tiles as being flagged and having a amplitude of 0. Default=[]",
    )

    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/mwa_utils",
        help="Dir where MWA FEE models will be saved. Default=./embers_out/mwa_utils",
    )

    _args = _parser.parse_args()
    _nside = _args.nside
    _pointings = [int(item) for item in _args.pointings.split(",")]

    _flags = [int(item) for item in _args.flags.split(",")]

    if _flags == [0]:
        _flags = []

    print(_flags)
    _out_dir = _args.out_dir

    print(f"MWA_FEE maps saved to: {_out_dir}")
    mwa_fee_model(_out_dir, _nside, _pointings, _flags)
