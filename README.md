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

&nbsp;

### Satellite Ephemeris  

72 satellites have been identified, which actively emit between 137.150 and 138.550 MHz. These satellites are identified by their unique NORAD Catalogue Number, which are hard-coded into `sat_ids.py`. The satellites of interest are ORBCOMM communication Satellites, NOAA and Meteor weather satellites. 

In this section of code we determine the trajectories of our satellites of interest.

`cd ./code/sat-ephemeris`

#### 1. Download Two Line Elements [ TLEs ]
Two Line Elements contain information about the orbital elements of Earth-orbiting objects. They can be used to determine the ephemeris of satellites. Relevant TLEs are downloaded from [Space-Track.org](https://www.space-track.org/) using the [spacetrack](https://pypi.org/project/spacetrack/) Python module. An account on space-track.org is required to download TLEs. 


`python download-tle.py --help` gives all options.  

`python download-tle.py --start_date=YYYY-MM-DD --stop_date=YYYY-MM-DD`   

#### 2. Satellite Ephemeris

A python astronomy package, [Skyfield](https://rhodesmill.org/skyfield/), is used to make sense of the TLEs. Every pair of lines (tle elements) in the TLE file corresponds to a particular epoch - a time at which the satellite's possition can be determined most accurately, before or after which it goes rapidly out of date. [Space-Track.org](https://www.space-track.org/) releases new TLEs for satellites at least once a day. `sat_ephemeris.py` determined the epochs of each pair of tle elements in a `TLE.txt` file, and used the elements which are most accurate. 

The MWA Satellite experiment is situated in remote Western Australia, at the Murchison Radio-astronomy Observatory. The latitude, longitude and elevation of the site are hard-coded into `sat_ephemeris.py`. Skyfield uses this in conjunction with the TLE files to determine the trajectories of satellite. `sat_ephemeris.py` computes the altitude and azimuth [alt/az] of the satellite at a particular time cadence. By default this is set to 20s, and this step is computationally expensive and can dramatically slow things down. 

```
python sat_ephemeris.py --help

python sat_ephemeris.py --sat=XXXXX --cadence=20 
```
This produces a rise time `t_rise`, set time `t_set`, in gps seconds, for each satellite pass. It also produces a corresponding list of altitudes and azimuths for each pass. These are saved to a json file, to be used later.


#### 3. Plot Ephemeris

Using the json files created in the previous section, we can now plot satellite passes. This will give us an idea of sky coverage, and how much integration is required to build up a complete beam map.

```
python plot_ephemeris.py --help

python plot_ephemeris.py --sat=XXXXX
```
&nbsp;

<p float="left">
  <img src="./documentation/21576.png" width="400" />
  <img src="./documentation/23545.png" width="400" /> 
</p>


<p float="left">
  <img src="./documentation/25981.png" width="400" />
  <img src="./documentation/25982.png" width="400" /> 
</p>


#### 4. Chronological Ephemeris

In later stages of this pipeline, we will need satellite ephemeris in chronological order. `chrono_ephem.py` collates data in all the json files produced by `sat_ephemeris.py`. Satellite passes are sorted by rise time, and saved to a new file `ultimate_ephem_list.json`.

```
python chrono_ephem.py --help

python chrono_ephem.py
```

### Beam Pointings

The MWA Telescope points its beam in different directions by giving each of its 16 dipoles a different delay. These pointings are quantised, with 181 discrete pointings available, each of which is given a unique 'Grid-pointing' number.  

It is important to know where the telescope is pointing at all times, because we will be making individual beam maps for each pointing. This metadata is available at [http://ws.mwatelescope.org](http://ws.mwatelescope.org/metadata/find). Metadata form the MWA telescope can be downloaded by following these steps.

```
cd ./code/beam-pointings/
```

let us download date between `2019-09-11 00:00:00` and `2019-11-02 00:00:00`. 
First, we need to convert these dates to gps-seconds.

```
python date-gps.py --date='2019-09-11 00:00:00'
gps: 1252195218

python date-gps.py --date='2019-11-02 00:00:00'
gps: 1256688018
```

Next, open [http://ws.mwatelescope.org](http://ws.mwatelescope.org/metadata/find). Enter the start and stop time. Also change the page size to 200. Leave all the other fields un-checked, as shown below. Finally, click `search`.



<p float="left">
  <img src="./documentation/metadata-1.jpg" width="1200" />
</p>

The important thing to note is the number of pages that are returned by the search. In this particular case it is 74.

<p float="left">
  <img src="./documentation/metadata-2.jpg" width="700" />
</p>

Using this information we can download the required metadata using `download-pointings.sh`. The site limits queries to approximately 200 per minute, which is why these steps were necessary.

```
source download_pointings.sh 74 1252195218 1256688018
```

The results of each page are downloaded to a json file in the `./../../ouputs/beam-pointings/` directory. The download script sleeps for 66s between each download, so as not to overload the MWA data servers.




