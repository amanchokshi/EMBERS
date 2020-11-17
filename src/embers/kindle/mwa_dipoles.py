"""
MWA Dipoles
-----------

Download MWA metafits files and determine flagged dipoles in the tiles

"""

import argparse

from embers.mwa_utils.mwa_dipoles import mwa_flagged_dipoles


def main():
    """
    Check MWA antenna dipole flagging using the :func:`~embers.mwa_utils.mwa_dipoles.mwa_flagged_dipoles` function.

    .. code-block:: console

        $ mwa_dipoles --help

    """

    _parser = argparse.ArgumentParser(
        description="""
        Download MWA metafits files and determine flagged dipoles in the tiles
        """
    )

    _parser.add_argument(
        "--num_files",
        metavar="\b",
        default=14,
        type=int,
        help="Number of metafits files to download Default: 20",
    )
    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/mwa_utils",
        help="Dir where MWA metadata will be saved. Default=./embers_out/mwa_utils",
    )

    _args = _parser.parse_args()
    _num_files = _args.num_files
    _out_dir = _args.out_dir

    mwa_flagged_dipoles(_num_files, _out_dir)
    print(f"MWA dipole flagging data saved to {_out_dir}")
