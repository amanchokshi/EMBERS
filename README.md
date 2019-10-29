## MWA Satellites
Code to measure the beam of MWA Tiles using satellites.

#### TLE
Two Line Elements contain information which can be used to determine the ephemeris of satellites at any given time. `download_TLE.py` contains a hardcoded dictionary of satellites of interest which emit between 137.150 and 138.550 MHz. TLEs are downloaded from [Space-Track.org](https://www.space-track.org/) using the [spacetrack](https://pypi.org/project/spacetrack/) Python module. An account on space-track.org is required to download TLEs. 

```
python download_TLE.py --start_date=YYYY-MM-DD --stop_date=YYYY-MM-DD
```

