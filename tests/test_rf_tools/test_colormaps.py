from pathlib import Path

from embers.rf_tools.colormaps import jade, plt_colormaps, spectral, waves_2d


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


def test_plt_colormaps():
    spec, spec_r = spectral()
    ja, ja_r = jade()
    plt_colormaps(spec, spec_r, ja, ja_r, ".")
    png = Path("colormaps.png")
    assert png.is_file() is True
    if png.is_file() is True:
        png.unlink()
