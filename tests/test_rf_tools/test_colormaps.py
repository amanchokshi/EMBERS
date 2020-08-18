from os import path
from pathlib import Path

from embers.rf_tools.colormaps import jade, spectral, waves_2d

# Save the path to this directory
dirpath = path.dirname(__file__)

# Obtain path to directory with test_data
test_data = path.abspath(path.join(dirpath, "../data"))


def test_spectral():
    spec, _ = spectral()
    assert type(spec).__name__ == "ListedColormap"


def test_spectral_r():
    _, spec_r = spectral()
    assert type(spec_r).__name__ == "ListedColormap"


def test_jade():
    ja, _ = jade()
    assert type(ja).__name__ == "ListedColormap"


def test_jade_r():
    _, ja_r = jade()
    assert type(ja_r).__name__ == "ListedColormap"


def test_waves_2d():
    wave = waves_2d()
    assert wave.shape == (512, 512)

