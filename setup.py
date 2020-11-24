"""Setup script for EMBERS package."""

from codecs import open
from pathlib import Path

from setuptools import find_packages, setup

# Get long description from README
with open("README.rst", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Get the requirements list
with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()


def data_files(data_dir, data_types):
    """Recursively return a list of data files."""
    data_dir = Path(data_dir)
    data_files = []
    for d in data_types:
        files = [p.relative_to(Path("embers")) for p in data_dir.rglob(d)]
        data_files.extend(files)
    return [p.as_posix() for p in data_files]


setup(
    name="embers",
    version="1.0.0",
    license="MIT",
    author="Aman Chokshi",
    author_email="achokshi@student.unimelb.edu.au",
    description="Experimental Measurement of BEam Response with Satellites",
    long_description=long_description,
    url="http://embers.readthedocs.io",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: MacOS",
        "Operating System :: Unix",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    entry_points={
        "console_scripts": [
            "colormaps=embers.kindle.colormaps:main",
            "waterfall_single=embers.kindle.waterfall_single:main",
            "waterfall_batch=embers.kindle.waterfall_batch:main",
            "align_single=embers.kindle.align_single:main",
            "align_batch=embers.kindle.align_batch:main",
            "download_tle=embers.kindle.download_tle:main",
            "ephem_single=embers.kindle.ephem_single:main",
            "ephem_batch=embers.kindle.ephem_batch:main",
            "ephem_chrono=embers.kindle.ephem_chrono:main",
            "sat_channels=embers.kindle.sat_channels:main",
            "mwa_pointings=embers.kindle.mwa_pointings:main",
            "mwa_dipoles=embers.kindle.mwa_dipoles:main",
            "mwa_fee=embers.kindle.mwa_fee:main",
            "ref_models=embers.kindle.ref_models:main",
            "rfe_calibration=embers.kindle.rfe_calibration:main",
            "tile_maps=embers.kindle.tile_maps:main",
            "null_test=embers.kindle.null_test:main",
            "compare_beams=embers.kindle.compare_beams:main",
        ],
    },
    keywords=("embers radio astronomy satellites beam measurement"),
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    include_package_data=True,
    package_data={"": data_files("embers/kindle/data", ["*.txt", "*.ffe"])},
    install_requires=requirements,
    python_requires=">=3.6, <4",
    zip_safe=False,
)
