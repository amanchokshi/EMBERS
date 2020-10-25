
===============
 API Reference
===============

Quick links to the sections below:

.. contents:: :local:

RF Tools
========

Contains :mod:`~embers.rf_tools.rf_data`, :mod:`~embers.rf_tools.align_data`, :mod:`~embers.rf_tools.colormaps` modules

.. currentmodule:: embers.rf_tools

.. autosummary::

    rf_data.read_data
    rf_data.tile_names
    rf_data.tile_pairs
    rf_data.time_tree
    rf_data.plt_waterfall
    rf_data.single_waterfall
    rf_data.batch_waterfall
    rf_data.waterfall_batch
    align_data.savgol_interp
    align_data.plot_savgol_interp
    align_data.save_aligned
    align_data.align_batch
    colormaps.spectral
    colormaps.jade
    colormaps.waves_2d
    colormaps.plt_colormaps


Sat Utils
=========

.. currentmodule:: embers.sat_utils

.. autosummary::

    sat_list.norad_ids
    sat_list.download_tle
    sat_ephemeris.load_tle
    sat_ephemeris.epoch_ranges
    sat_ephemeris.epoch_time_array
    sat_ephemeris.sat_pass
    sat_ephemeris.ephem_data
    sat_ephemeris.sat_plot
    sat_ephemeris.save_ephem
    sat_ephemeris.ephem_batch
    chrono_ephem.obs_times
    chrono_ephem.interp_ephem
    chrono_ephem.write_json
    chrono_ephem.save_chrono_ephem
    sat_channels.read_aligned
    sat_channels.noise_floor
    sat_channels.time_filter
    sat_channels.plt_window_chans
    sat_channels.plt_channel
    sat_channels.plt_sats
    sat_channels.good_chans
    sat_channels.window_chan_map
    sat_channels.batch_window_map


MWA Utils
=========

.. currentmodule:: embers.mwa_utils

.. autosummary::

    mwa_pointings.download_meta
    mwa_pointings.clean_meta_json
    mwa_pointings.combine_pointings
    mwa_pointings.point_integration
    mwa_pointings.pointing_hist
    mwa_pointings.rf_obs_times
    mwa_pointings.obs_pointings
    mwa_pointings.tile_integration
    mwa_pointings.plt_hist_array
    mwa_pointings.mwa_point_meta
    mwa_dipoles.download_metafits
    mwa_dipoles.find_flags
    mwa_dipoles.mwa_flagged_dipoles
    mwa_fee.local_beam
    mwa_fee.mwa_fee_model


Tile Maps
=========

   .. currentmodule:: embers.tile_maps

.. autosummary::

    beam_utils.rotate_map
    beam_utils.healpix_cardinal_indices
    beam_utils.healpix_cardinal_slices
    beam_utils.nan_mad
    beam_utils.map_slices
    beam_utils.poly_fit
    beam_utils.chisq_fit_gain
    beam_utils.chisq_fit_test
    beam_utils.plt_slice
    beam_utils.plot_healpix
    ref_fee_healpix.create_model
    ref_fee_healpix.ref_healpix_save
    tile_maps.check_pointing
    tile_maps.plt_channel
    tile_maps.plt_fee_fit
    tile_maps.rf_apply_thresholds
    tile_maps.rfe_calibration
    tile_maps.rfe_collate_cali
    tile_maps.rfe_batch_cali
    tile_maps.project_tile_healpix
    tile_maps.mwa_clean_maps
    tile_maps.plt_sat_maps
    tile_maps.plt_clean_maps
    tile_maps.tile_maps_batch
    null_test.good_ref_maps
    null_test.plt_null_test
    null_test.null_test
    compare_beams.beam_slice
    compare_beams.batch_compare_beam


Kindle
======
