# MWA Satellites
Code to measure the beam shapes of MWA Tiles using satellites.  

### Conda Environment

A conda virtual enviroment is used to install the correct version of python and all required packages in an isolated enviroment, ensuring that our code performs as expecteted. Either [Anaconda](https://www.anaconda.com/distribution/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) are required before proceeding further. All the required packages are in `env.yml`, which can be installed using:
```
source setup.sh
```
This creates a conda enviroment, called sat-env, for out mwa-satellites repository. This virtual conda enviroment can be activated and deactivated with:

```
conda activate sat-env

conda deactivate
```

When running code in this repository, ensure that `sat-env`is active.


### Satellite Ephemeris  

72 satellites have been identified, which actively emit between 137.150 and 138.550 MHz. These satellites are identified by their unique NORAD Catalogue Number, which are hard-coded into `sat_ids.py`. The satellites of interest are ORBCOMM communication Satellites, NOAA and Meteor weather satellites. 

In this section of code we determine the trajectories of our satellites of interest.

`cd ./code/sat-ephemeris`

#### 1. Download Two Line Elements [ TLEs ]
Two Line Elements contain information about the orbital elements of Earth-orbiting objects. They can be used to determine the ephemeris of satellites. Relevant TLEs are downloaded from [Space-Track.org](https://www.space-track.org/) using the [spacetrack](https://pypi.org/project/spacetrack/) Python module. An account on space-track.org is required to download TLEs. 


`python download-tle.py --help` gives all options.  

`python download-tle.py --start_date=YYYY-MM-DD --stop_date=YYYY-MM-DD`   

#### 2. Satellite Ephemeris

A python astronomy package, [Skyfield](https://rhodesmill.org/skyfield/), is used to make sense of the TLEs. Every pair of lines (tle elements) in the TLE file corresponds to a particular epoch - a time at which the satellite's possition can be determined most accurately, before or after which they go rapidly out of date. [Space-Track.org](https://www.space-track.org/) releases new TLEs for satellites at least once a day. `sat_ephemeris.py` determined the epochs of each pair of tle elements in the `TLE.txt` and uses this to construct a time array. Intervals between elements of this time array correspond to a time range when a particular pair of tle elements is most accurate. Skyfield deals with time most comfortably in its julian date format. [Astropy](https://www.astropy.org/) is used to convert this to gps seconds. 








