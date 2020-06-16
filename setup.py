# Package imports
from setuptools import setup, find_namespace_packages
import re

# Get long description from README
with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()

# Get the requirements list
with open('requirements.txt', 'r') as f:
    requirements = f.read().splitlines()
    requirements = requirements[:-1]

# Read the __version__.py file
with open('embers/__version__.py', 'r') as f:
    vf = f.read()

# Obtain version from read-in __version__.py file
version = re.search(r"^_*version_* = ['\"]([^'\"]*)['\"]", vf, re.M).group(1)

setup(
    name="embers",
    version="0.0.1",
    license='MIT',
    author="Aman Chokshi",
    author_email="achokshi@student.unimelb.edu.au",
    description="Experimental Measurement of BEam Response with Satellites",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/amanchokshi/EMBERS",
    packages=find_namespace_packages(include=["embers.*"]),
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
    keywords=("embers radio astronomy satellites beam measurement"),
    python_requires='>=3.6, <4',
    install_requires=requirements,
    dependency_links=['https://github.com/MWATelescope/mwa_pb@master#egg=mwa_pb'],
    zip_safe=False,
)