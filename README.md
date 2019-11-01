## MWA Satellites
Code to measure the beam of MWA Tiles using satellites.

#### 1. Two Line Elements [ TLEs ]
Two Line Elements contain information the orbital elements of Earth-orbiting objects. TLEs can be used to determine the ephemeris of satellites. For this experimet, 72 satellites have been identified, which actively emit between 137.150 and 138.550 MHz. These satellites are identified by their unique NORAD Catalogue Number, which are hard-coded into `download_TLE.py`. The satellites of interest are ORBCOMM communication Satellites, NOAA and Meteor weather satellites. Relevant TLEs are downloaded from [Space-Track.org](https://www.space-track.org/) using the [spacetrack](https://pypi.org/project/spacetrack/) Python module. An account on space-track.org is required to download TLEs.

```
python download_TLE.py --start_date=YYYY-MM-DD --stop_date=YYYY-MM-DD
```   

#### 2. Satellite Ephemeris



