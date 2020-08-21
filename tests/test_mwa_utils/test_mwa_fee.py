import json
import shutil
from os import path
from pathlib import Path

import healpy as hp
import numpy as np
from embers.mwa_utils.mwa_fee import local_beam, mwa_fee_model

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))

nside = 32
npix = hp.nside2npix(nside)
npix = hp.nside2npix(nside)
beam_response_XX = np.zeros(npix)
beam_response_YY = np.zeros(npix)
above_horizon = range(int(npix / 2))
beam_zas, beam_azs = hp.pix2ang(nside, above_horizon)

response = local_beam([list(beam_zas)], [list(beam_azs)], freq=137e6, power=False)
print(response[0][0][0])


def test_local_beam():
    response = local_beam([list(beam_zas)], [list(beam_azs)], freq=137e6, interp=False,)
    assert response[0][0][0] == 0.9911146621736204


def test_local_beam_interp():
    response = local_beam([list(beam_zas)], [list(beam_azs)], freq=137e6, interp=True,)
    assert response[0][0][0] == 0.9910795273080615


def test_local_beam_jones():
    response = local_beam([list(beam_zas)], [list(beam_azs)], freq=137e6, jones=True,)
    assert response[0][0][0][0] == (0.7036201076403841 - 0.011709834111186911j)


def test_local_beam_power_false():
    response = local_beam([list(beam_zas)], [list(beam_azs)], freq=137e6, power=False,)
    assert response[0][0][0] == 0.9955297721856747
