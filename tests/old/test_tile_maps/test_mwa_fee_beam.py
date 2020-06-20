import sys
import unittest
import numpy as np
import healpy as hp
from pathlib import Path

sys.path.append("../../code/tile_maps")
import mwa_fee_beam as mfb
from mwa_pb.mwa_sweet_spots import all_grid_points


class Test_mwa_fee_beam(unittest.TestCase):

    nside = 32
    npix = hp.nside2npix(nside)

    # healpix indices above horizon
    # convert to zenith angle and azimuth
    above_horizon = range(int(npix / 2))
    beam_zas, beam_azs = hp.pix2ang(nside, above_horizon)
    delay_p = np.array([all_grid_points[int("0")][-1], all_grid_points[int("0")][-1]])
    amps = np.ones((2, 16))

    def test_create_model(self):
        response = mfb.local_beam(
            [list(self.beam_zas)],
            [list(self.beam_azs)],
            freq=137e6,
            delays=self.delay_p,
            zenithnorm=True,
            power=True,
            interp=False,
            amps=self.amps,
        )

        response_XX = response[0][0]
        response_YY = response[1][0]

        self.assertEqual(response_XX.shape, response_YY.shape)
        self.assertEqual(response_XX[0], 0.9911146621736204)
        self.assertEqual(response_YY[0], 0.9911152950440345)


if __name__ == "__main__":
    unittest.main()
