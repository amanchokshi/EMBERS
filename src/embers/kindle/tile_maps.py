"""
Tile Maps
=========

Batch process satellite RF data to create MWA beam maps and intermediate plots.
Outputs saved to ``./embers_out/tile_maps/tile_maps``

"""

import argparse

from embers.tile_maps.tile_maps import tile_maps_batch


def main():
    """
    Create tile beam maps using the :func:`~embers.tile_maps.tile_maps.tile_maps_batch`.

    .. code-block:: console

        $ tile_maps --help

    """

    _parser = argparse.ArgumentParser(
        description="""
        Batch process satellite RF data to create clean beam maps and all intermediate data products.
        """
    )

    _parser.add_argument(
        "--start_date",
        metavar="\b",
        required=True,
        help="start date in YYYY-MM-DD format",
    )

    _parser.add_argument(
        "--stop_date",
        metavar="\b",
        required=True,
        help="stop date in YYYY-MM-DD format",
    )

    _parser.add_argument(
        "--sat_thresh",
        metavar="\b",
        default=1,
        type=int,
        help="σ threshold to detect sats in the computation of rf data noise_floor. Default: 1",
    )

    _parser.add_argument(
        "--noi_thresh",
        metavar="\b",
        default=3,
        type=int,
        help="noise threshold: multiples of mad. default: 3",
    )

    _parser.add_argument(
        "--pow_thresh",
        metavar="\b",
        default=5,
        type=int,
        help="Peak power which must be exceeded for satellite pass to be considered. default: 5dB",
    )

    _parser.add_argument(
        "--ref_model",
        metavar="\b",
        default="embers_out/tile_maps/ref_models/ref_dipole_models.npz",
        help="Path to reference feko model. Default: embers_out/tile_maps/ref_models/ref_dipole_models.npz",
    )

    _parser.add_argument(
        "--fee_map",
        metavar="\b",
        default="embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
        help="Path to MWA FEE model. Default: embers_out/mwa_utils/mwa_fee/mwa_fee_beam.npz",
    )

    _parser.add_argument(
        "--rfe_cali",
        metavar="\b",
        default="embers_out/tile_maps/rfe_calibration/rfe_gain_fit.npy",
        help="Path to RF Explorer gain calibration solution: embers_out/tile_maps/rfe_calibration/rfe_gain_fit.npy",
    )

    _parser.add_argument(
        "--nside", metavar="\b", default=32, type=int, help="Healpix nside. Default: 32"
    )

    _parser.add_argument(
        "--obs_point_json",
        metavar="\b",
        default="embers_out/mwa_utils/obs_pointings.json",
        help="Path to obs_pointings.json. Default: embers_out/mwa_utils/obs_pointings.json",
    )

    _parser.add_argument(
        "--align_dir",
        metavar="\b",
        default="embers_out/rf_tools/align_data",
        help="Directory where aligned rf data lives. Default: embers_out/rf_tools/align_data",
    )

    _parser.add_argument(
        "--chrono_dir",
        metavar="\b",
        default="embers_out/sat_utils/ephem_chrono",
        help="Directory where Chronological satellite ephemeris data lives. Default: embers_out/sat_utils/ephem_chrono",
    )

    _parser.add_argument(
        "--chan_map_dir",
        metavar="\b",
        default="embers_out/sat_utils/sat_channels/window_maps",
        help="Directory where satellite frequency channel maps are saved. Default: embers_out/sat_utils/sat_channels/window_maps",
    )

    _parser.add_argument(
        "--out_dir",
        metavar="\b",
        default="./embers_out/tile_maps/tile_maps",
        help="Directory where RF Explorer calibration data will be saved. Default=./embers_out/tile_maps/tile_maps",
    )

    _parser.add_argument(
        "--plots",
        metavar="\b",
        default="False",
        help="If True, create a zillion diagnostic plots at the project_tile_healpix stage. Default: False",
    )

    _parser.add_argument(
        "--rfe_cali_bool",
        metavar="\b",
        default="True",
        help="If True, apply RFE calibration. Default: True",
    )

    _parser.add_argument(
        "--max_cores",
        metavar="\b",
        type=int,
        help="Maximum number of cores to be used by this script. By default all core available cores are used",
    )

    _args = _parser.parse_args()
    _start_date = _args.start_date
    _stop_date = _args.stop_date
    _sat_thresh = _args.sat_thresh
    _noi_thresh = _args.noi_thresh
    _pow_thresh = _args.pow_thresh
    _ref_model = _args.ref_model
    _fee_map = _args.fee_map
    _rfe_cali = _args.rfe_cali
    _nside = _args.nside
    _obs_point_json = _args.obs_point_json
    _align_dir = _args.align_dir
    _chrono_dir = _args.chrono_dir
    _chan_map_dir = _args.chan_map_dir
    _out_dir = _args.out_dir
    _plots = _args.plots
    _rfe_cali_bool = _args.rfe_cali_bool
    _max_cores = _args.max_cores

    if _plots == "False":
        _plots is False
    else:
        _plots = True

    if _rfe_cali_bool == "True":
        _rfe_cali_bool is True
    else:
        _rfe_cali_bool is False

    tile_maps_batch(
        _start_date,
        _stop_date,
        _sat_thresh,
        _noi_thresh,
        _pow_thresh,
        _ref_model,
        _fee_map,
        _rfe_cali,
        _nside,
        _obs_point_json,
        _align_dir,
        _chrono_dir,
        _chan_map_dir,
        _out_dir,
        _plots,
        _rfe_cali_bool,
        max_cores=_max_cores,
    )

    print(f"MWA tile map files saved to: {_out_dir}")
