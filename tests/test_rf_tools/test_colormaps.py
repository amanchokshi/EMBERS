import pytest
from embers.rf_tools.colormaps import spectral, jade


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
