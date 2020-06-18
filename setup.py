# Package imports
from setuptools import setup, find_packages
from pathlib import Path
import re

# Get long description from README
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Get the requirements list
with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

# Recursive data files in embers.kindle
kindle_data = Path(f"{__file__}/embers/kindle/data")
files = [str(p.relative_to(kindle_data)) for p in kindle_data.rglob("*.txt")]

# Obtain name, version, author from __init__.py file
from embers import __name__, __version__, __author__

name = __name__
version = __version__
author = __author__


setup(
    name=name,
    version=version,
    license="MIT",
    author=author,
    author_email="achokshi@student.unimelb.edu.au",
    description="Experimental Measurement of BEam Response with Satellites",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/amanchokshi/EMBERS",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
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
    entry_points={"console_scripts": ["colormaps=embers.kindle.colormaps:main"],},
    keywords=("embers radio astronomy satellites beam measurement"),
    python_requires=">=3.6, <4",
    install_requires=requirements,
    include_package_data=True,
    zip_safe=False,
)
