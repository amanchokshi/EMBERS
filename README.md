# MWA Satellites
**E**xperimental **M**easurement of **BE**am **R**esponses with **S**atellites  

### Installation

The *EMBERS* package may be installed by cloning the repository. Using a virtual environment is suggested to prevent conflicting dependancies.

```
git clone https://github.com/amanchokshi/EMBERS.git
cd EMBERS
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```


&nbsp;
#### Directory Structure

```
.
├── code                    # Code lives here
├── docs                    # Documentation files
├── data                    # Sample data and beam maps
├── outputs                 # All outputs saved here by default
├── setup.sh                # Creates conda virtual env with dependancies
├── env.yml                 # Conda environment file
├── LICENSE                 
└── README.md
```

&nbsp;
### Decode RF Explorer Data

The RF Explorers saves a single line of data, per time step, to a text file. This contains a timestamp, in UNIX time, and a binary string of power data, corresponding to power in dB per frequency channel. `rf_data.py` has two functions. The first, decodes the rf data files and returns a list of times, and a 2D array of powers. The second function used the time and power arrays to create a waterfall plot. This is an extremely useful diagnostic tool, displaying satellite passes at bright vertical streaks. 

```
cd ./code/decode_rf_data/

python batch_waterfall.py --help
python batch_waterfall.py --data_dir=../../data/rf_data/ --start_date=2019-10-10 --stop_date=2019-10-10
```

Can be used to create a batch of waterfall plots from data in a time interval. The above code creates waterfall plots from the two sample data files provided in the `data/rf_data` directory, saving the results in `outputs/decode_rf_data`. On the left is the waterfall plot of the reference antenna, with an obvious broadband weather satellite pass visible at the bottom right. The narrower vertical lines are ORBCOMM satellites. The plot on the right shows the same satellites as seen by an MWA tile. The nulls in the primary beam are clearly visible. 

<p float="left">
  <img src="./imgs/rf0XX_2019-10-10-02:30.png" width="49%" />
  <img src="./imgs/S06XX_2019-10-10-02:30.png" width="49%" /> 
</p>

A `logs.txt` file is saved to the output directory. This shows what plots were saved by `batch_waterfall.py`, and also which data files were missing.

&nbsp;
### Align Data

The RF Explorers do not record data at exactly the same frequency and do not start recording at exactly the same time. In fact, the older models record at approximately 6 Hz, while the newer ones are capable of a sampling rate of nearly 9 Hz. This discrepency in sampling rates makes it difficult to compare any two data samples. This issue is overcome by smoothing the data, along the time axis, with a Savitzky-Golay filter. Interpolating the smoothed data and resampling it at a constant frequency [ 0.5 Hz ] gives us a easier data set to work with. The starting times of the power arrays are also synchronized.

```
cd ../align_data

python align_data.py
```

<p float="left">
  <img src="./imgs/align_data.png" width="110%" /> 
</p>

This illustrates how `align_data.py` smoothes the noisy power data. The image displayed above represents only one frequency channel  [59] from the 2D power array. This channel was chosen to best show the spectral features of satellite passes as seen by the reference and mwa tiles. 

In practince, savgol smoothing and interpolation are applied to the 2D power arrays, along the time axis.

Aligned pairs of data sets are saved using `batch_align.py`. Using a start date, and stop date, it saves aligned data sets as .npz files, for all raw data found in the time interval.

```
python batch_align.py --help
```

A `logs.txt` file is saved to the output directory. This shows which data sets were saved by `batch_align.py`, and also which data files were missing.


&nbsp;
### Satellite Ephemeris  

72 satellites have been identified, which actively emit between 137.150 and 138.550 MHz. These satellites are identified by their unique NORAD Catalogue Number, which are hard-coded into `sat_ids.py`. The satellites of interest are ORBCOMM communication Satellites, NOAA and Meteor weather satellites. 

In this section of code we determine the trajectories of our satellites of interest.

`cd ../sat_ephemeris`

#### 1. Download Two Line Elements [ TLEs ]
Two Line Elements contain information about the orbital elements of Earth-orbiting objects. They can be used to determine the ephemeris of satellites. Relevant TLEs are downloaded from [Space-Track.org](https://www.space-track.org/) using the [spacetrack](https://pypi.org/project/spacetrack/) Python module. An account on space-track.org is required to download TLEs. 


`python download_tle.py --help` gives all options.  

`python download_tle.py --start_date=YYYY-MM-DD --stop_date=YYYY-MM-DD`   

#### 2. Satellite Ephemeris

A python astronomy package, [Skyfield](https://rhodesmill.org/skyfield/), is used to make sense of the TLEs. Every pair of lines (tle elements) in the TLE file corresponds to a particular epoch - a time at which the satellite's possition can be determined most accurately, before or after which it goes rapidly out of date. [Space-Track.org](https://www.space-track.org/) releases new TLEs for satellites at least once a day. `sat_ephemeris.py` determined the epochs of each pair of tle elements in a `TLE.txt` file, and used the elements which are most accurate. 

The MWA Satellite experiment is situated in remote Western Australia, at the Murchison Radio-astronomy Observatory. The latitude, longitude and elevation of the site are hard-coded into `sat_ephemeris.py`. Skyfield uses this in conjunction with the TLE files to determine the trajectories of satellite. `sat_ephemeris.py` computes the altitude and azimuth [alt/az] of the satellite at a particular time cadence. By default this is set to 2 seconds, and this step is computationally expensive and can dramatically slow things down. 

```
python sat_ephemeris.py --help

python sat_ephemeris.py --sat=XXXXX --cadence=20 
```
This produces a rise time `t_rise`, set time `t_set`, in gps seconds, for each satellite pass. It also produces a corresponding list of altitudes and azimuths for each pass. These are saved to a json file, to be used later. The following can be used to create ephemeris json files for all satellite tle files.

```
python batch_sat_ephem.py --help
```

A log file is also saved in the `./../../output/sat_ephemeris/ephem_json/sat_ephem_logs.txt` directory, showing if any errors occured.


&nbsp;
#### 3. Plot Ephemeris

Using the json files created in the previous section, we can now plot satellite passes. This will give us an idea of sky coverage, and how much integration is required to build up a complete beam map.

```
python plot_ephemeris.py --help

python plot_ephemeris.py
```
This will plot all the json ephemeris files that were generated in the previous section, and will save them to `../../outputs/sat_ephemeris/ephem_plots/`

&nbsp;

<p float="left">
  <img src="./imgs/sat_ephem.jpg" width="100%" />
</p>


#### 4. Chronological Ephemeris

In later stages of this analysis we will need satellite ephemeris in chronological order. `chrono_ephem.py` collates data in all the json files produced by `sat_ephemeris.py`. As a first step, the gps timestamps used by `skyfield`, are converted to unix timestamps, to match the format used in the rf explorer data. Next, the satellite alt/az ephemeris are interpolated, to match the time cadence of `align_data.py`. It is crutial to remember to set `--interp_freq` to be the same as that used in `align_data.py`, or just use the default of 0.5 Hz. For each 30 min observation, a corresponding sat ephem file is created in `./../../outputs/sat_ephemeris/chrono_json/` which contains all the satellite passes within the given time window. 

```
python chrono_ephem.py --help
```


&nbsp;
### Satellite Channels

As access to the ORBCOMM interface box was not available, the channels in which each satellite transimits can be determined with a careful analysis of the RF data and satellite ephemeris.

We use reference data to detect satellite channels because it has the best SNR. Pairing a reference RF data file, with it's corresponding chrono_ephem json file gives us the satellite expected within each 30 minute observation. Looping over the satellites in the chrono_ephem files, we identify the temporal region of the rf data where we expect to see its signal. Let us call this the *window*. We now use a series of thresholding criteria to help identify the most probable channel.

#### Channel Filtering 

The following thresholds were used to identify the correct channel:

* Power threshold
* Noise threshold
* Window Occupancy


The peak signal in our data is `~ 30 dB` above the noise floor. We define an power threshold at `15 dB` and require that the maximum signal in a channel must exceed this value, if it contains a satellite.

Our observation windows are sparesely populated with satellite signal. We classify any channel with a maximum signal below *1σ* to be *noisy data*. We find the *median [μ_noise]* of the *noisy data* and it *Median Absolute Deviatio [σ_noise]*. Using these we construct a noise threshold of *μ_noise + N x σ_noise*. By default N = 3. We now use this noise threshold to more rigorously choose channels with satellites present.

We also demand that the signal at the edges of the *window* be below the noise threshold. This ensures that the satellite begins and ends within the expected region. 

A fractional occupancy of the window is computed by finding length of signal above the noise floor and dividing it by the window length. We require that *0.8 ≤ occupancy < 1.0*. By requiring that occupancy is less than `1.0`, we ensure that satellite passes longer than the window are not classified as potential channels.

The satellites periodically shift frequency channel over the 5 month period of our experiment, we make a *window map* for every half hour observation using `window_chan_map.py`. Each map identified the channels of all satellites in the 30 minute observation using the criteria described above - noise, power and occupation thresholds. if the `--plots=True` argparse option is used, a waterfall plot of each observation are generated, with the time window of the satellite and possible channels highlighted. The most probable channel is highlighted in green, as seen below.

```
cd ../sat_channels

python window_chan_map.py --help
```


<p float="left">
  <img src="./imgs/2019-10-10-02:30_40086_[4, 8]_waterfall.png" width="49%" />
  <img src="./imgs/2019-10-10-02:30_44387_[56, 57, 58, 59, 60, 61, 62, 63]_waterfall.png" width="49%" />
</p>

<p float="left">
  <img src="./imgs/2019-10-10-02:30_40086_4_channel.png" width="49%" />
  <img src="./imgs/2019-10-10-02:30_44387_59_channel.png" width="49%" />
</p>


&nbsp;
### Reproject Reference Beam Models

FEKO simulations were run to generate FEE beam models for the reference antennas. A model was generated for the XX & YY polarization of the reference antennas and is saved in `/data/FEE_Reference_Models/`. The `.ffe` files generated by FEKO need to be converted to something that we can use. `project_beam2healpix.py` finds the power in the beam, interpolates and rotates it into a healpix projection. The FEE beam is scaled such that its power at zenith is 0 dB.

```
cd ../reproject_ref

python project_beam2healpix.py
```

This will create `ref_dipole_models.npz` and a diagnostic plot of the beams in the output directory `../../outputs/reproject_ref/reproject_dipole_models.png`.

<p float="left">
  <img src="./imgs/reproject_dipole_models.png" width="110%" /> 
</p>



&nbsp;
### Beam Pointings

The MWA Telescope points its beam in different directions by giving each of its 16 dipoles a different delay. These pointings are quantised, with 181 discrete pointings available, each of which is given a unique 'Grid-pointing' number.  

It is important to know where the telescope is pointing at all times, because we will be making individual beam maps for each pointing. This metadata is available at [http://ws.mwatelescope.org](http://ws.mwatelescope.org/metadata/find). Metadata form the MWA telescope can be downloaded by following these steps.

```
cd ../beam-pointings/
```

let us download date between `2019-09-11 00:00:00` and `2019-11-02 00:00:00`. 
First, we need to convert these dates to gps-seconds.

```
python date_gps.py --date='2019-09-11 00:00:00'
gps: 1252195218

python date_gps.py --date='2019-11-02 00:00:00'
gps: 1256688018
```

Next, open [http://ws.mwatelescope.org](http://ws.mwatelescope.org/metadata/find). Enter the start and stop time. Also change the page size to 200. Leave all the other fields un-checked, as shown below. Finally, click `search`.



<p float="left">
  <img src="./imgs/metadata-1.jpg" width="1200" />
</p>

The important thing to note is the number of pages that are returned by the search. In this particular case it is 74.

<p float="left">
  <img src="./imgs/metadata-2.jpg" width="700" />
</p>

Using this information we can download the required metadata using `download_pointings.sh`. The site limits queries to approximately 200 per minute, which is why these steps were necessary.

```
source download_pointings.sh 74 1252195218 1256688018
```

The results of each page are downloaded to a json file in the `./../../ouputs/beam_pointings/` directory. The download script sleeps for 66s between each download, so as not to overload the MWA data servers.


The next step is to collate the data in all these discrete json metadata files into and `ultimate_pointing_times.json` list. This is slightly tricky, as there are a days when the MWA is idle as a whole, but the MWA Satellite experiment is recording data. This means that the metadata we downloaded will not contain any info on these times, and we must infer the pointing of the telescope from the last pointing the receivers were in. 

```
python sort_pointings.py
```

This creates an `ultimate_pointing_times.json` file in the metadata directory. `pointing_list.py` has condensed observations with the same pointing into blocks of observation. Each entry in `ultimate_pointing_times.json` contains the `grid_pt` number, `start_gps`, `stop_gps` time and `obs_length` in seconds.

A good way to visualise the amount of data collected so far is to use `plot_pointings.py`. This computes total integration over each pointing, and converts it into hours, and saves a plot to `pointing_integration.png`. It plots all pointings with a tottal observation time above a threshold, in this case 200 hours. The plot below represent data between `2019-09-10` and `2020-03-17`. We see that pointings `0, 2, 4, 41` have more than 200 hours each, and we will make maps of each of them.

```

python plot_pointings.py
```

<p float="left">
  <img src="./imgs/pointing_integration.png" width="100%" />
</p>

We now need to classify the pointing of each 30 minute observation. 

```
python obs_pointings.py --help
```

This creates `obs_pointings.json` in the output directory, which contain three lists, for the three pointings which we are interested in. Each list contains the observations with the correct pointings. Next, using `tile_pointings.py` we create `tile_pointing_integration.png`, in the outputs folder, to examine how much integration each of the 32 tiles have. Because we had continuous failures of our equipment, some tiles have significanly less data than others.

```
python tile_pointings.py --help
```

<p float="left">
  <img src="./imgs/tiles_pointing_integration.png" width="100%" />
</p>

We also need to check whether any of the tiles have flagged dipoles. This information can be found in the metafits files saved by the mwa. We download these files with

```
python download_metafits.py
```

Extract data for the 14 tiles used in the experiment. [Hardcoded list]. Find Which of their dipoles are flagged and when. A value of `0`, indicates 0 dipoles flagged, while a value of `1-16` corresponds to the position of each of the 16 dipoles.

```
python dead_dipoles.py
```

<p float="left">
  <img src="./imgs/dead_dipoles.png" width="98%" />
</p>

We can see that dipole `9` of tile `HexS33Y` is flagged. We need to take this into account in later stages.


&nbsp;
### Tile Maps

In this section we create and analyse tile maps. We compare out tile maps with the latest FEE simulated beam models of MWA tiles. To do this, we need an external library [https://github.com/MWATelescope/mwa_pb.git](https://github.com/MWATelescope/mwa_pb.git), and the raw FEE Beam model data which is stored on cerberus [http://cerberus.mwa128t.org]. These were downloaded and installed while setting up the conda environment.

```
cd ../tile_maps
```

We can then use `mwa_fee_beam.py` to make simulated FEE beam maps at the 4 pointings of interest. This created plots and a `.npz` file with the four beam maps in the `/outputs/tile_analysis/FEE_maps directory`.

#### MWA FEE Beam Maps
<p float="left">
  <img src="./imgs/mwa_fee_beam_4.png" width="24%" />
  <img src="./imgs/mwa_fee_beam_0.png" width="24%" />
  <img src="./imgs/mwa_fee_beam_2.png" width="24%" />
  <img src="./imgs/mwa_fee_beam_41.png" width="24%" />
</p>

Finally, everything can come together.  We project the aligned satellite data onto healpix maps according to their ephemeris, using the channel maps to only project the good data. We use the power and noise threshold described in the satellite channel section. Signal which pass both thresholds are considered good, and their corresponding altitude and azimuth are determined from the ephem file. The RF Explorers distort input signals above -30dBm. We clip any such data to prevent inaccuracies in out final maps. Each satellite pass is git to the fee model using a single gain offset, by a least squares method. A goodness of fit chi-squared test is also used to reject outliers which may have slipped throught he channel map section. The resulting normalized passes are projected to a healpix map. The crutial step of dividing the tile data by the reference data, and then multyplying by the reference fee model is also done using `tile_maps.py` The resulting maps are stored in `.npz` files in the `./outputs/tile_maps/tile_maps_raw` directory. The data is stored as nested dictionaries.

```
cd ../tile_maps/

python tile_maps --help
```

The resulting `.npz` files have the following internal structure.

```
data
├── mwa_map                   
│   └── pointings            
│       └── satellites            
│           └── healpix maps            
├── ref_map                   
│   └── pointings            
│       └── satellites            
│           └── healpix maps            
├── tile_map                  
│   └── pointings            
│       └── satellites            
│           └── healpix maps            
└── time_map
    └── pointings            
        └── satellites            
            └── healpix maps
```

The highest level dictionary contains normalized mwa, reference, tile and time maps. Within each of these, there are dictionaries for each of the telescope pointings:`0, 2, 4, 41`. Within which there are dictionaries for each satellite norad ID, which contain a healpix map of data from one satellite, in one pointing. This structure map seem complicated, but is very useful for diagnostic purposes, and determining where errors in the final tile maps come from. The time maps contain the times of every data point added to the above maps.


Using `plot_sat_maps.py`, we can plot the sky coverage of each satellite. This is very useful in limiting the amount of "bad" data we add to our final maps. Some of the maps in `./outputs/tile_maps/sat_maps/` clearly show sparse coverage and weird power structure. Eyeballing these maps, we create a list of "good_sats", which we will use to make our good maps.

```
python plot_sat_maps.py
```

#### Satellite Maps
<p float="left">
  <img src="./imgs/23545_0_S07XX_rf0XX_passes.png" width="32%" />
  <img src="./imgs/28654_0_S07XX_rf0XX_passes.png" width="32%" />
  <img src="./imgs/41188_0_S07XX_rf0XX_passes.png" width="32%" />
</p>

These are three example sat maps. The first map has almost no data, so we don't use it in our final map. The other two show a lot of satellite passes / sky coverage, and we can almost make out the beam structure. The satellites will be included in out final "good maps". 

If a subset of satellite need to be used, they can be added to the `good_sats` list in `tile_maps_norm.py`. By default all satellites are used. We can now proceed to extracting and analysing all our "good" data from the `./outputs/tile_maps/tile_maps_raw`. `tile_maps_norm.py`, extracts the good satellite reference and tile data. This gives us the measured beam shape of the MWA tile, which are saved to `./outputs/tile_maps/tile_maps_norm/`.

```
python tile_maps_norm.py
```


We can now finally create some beam maps using `plot_tile_maps.py`, which caluclates the median of sat passes in each healpix pixel. The Maps are saved to `../../outputs/tile_maps/good_maps`, where good tile maps are created for each pointing, and for every tile. Error and count maps are also created at the same time.

```
python plot_tile_maps.py
```


#### Zenith Maps
<p float="left">
  <img src="./imgs/S35XX_rf0XX_0_good_map.png" width="32%" />
  <img src="./imgs/S35XX_rf0XX_0_good_map_counts.png" width="32%" />
  <img src="./imgs/S35XX_rf0XX_0_good_map_errors.png" width="32%" />
</p>

#### Sweet Point 2
<p float="left">
  <img src="./imgs/S35XX_rf0XX_2_good_map.png" width="32%" />
  <img src="./imgs/S35XX_rf0XX_2_good_map_counts.png" width="32%" />
  <img src="./imgs/S35XX_rf0XX_2_good_map_errors.png" width="32%" />
</p>

#### Sweet Point 4
<p float="left">
  <img src="./imgs/S35XX_rf0XX_4_good_map.png" width="32%" />
  <img src="./imgs/S35XX_rf0XX_4_good_map_counts.png" width="32%" />
  <img src="./imgs/S35XX_rf0XX_4_good_map_errors.png" width="32%" />
</p>

#### Sweet Point 41
<p float="left">
  <img src="./imgs/S35XX_rf0XX_41_good_map.png" width="32%" />
  <img src="./imgs/S35XX_rf0XX_41_good_map_counts.png" width="32%" />
  <img src="./imgs/S35XX_rf0XX_41_good_map_errors.png" width="32%" />
</p>


#### Null Test

We also perform a sanity check, to make sure that both our reference antennas have the same performance. We do this by fitting slices of the reference beam model to corresponding slices of out projected healpix data. We fit the residuals with a 3rd order polinomial to show what it has no structure. We then divide corresponding slices from `rf0 & rf1`


Perform the actual null test and create some interesting plots with

```
python null_test.py --help
```

We now compare corresponding EW and NS slices of both reference antennas. The four reference maps we created above are slices along the two cardinal axes. We compute the meadian value of each pixel, and estimate it's error using the Median Absolute Deviation. We rotate our FEE models by (π/4), to match the rotation of our data, and slice it along NS & EW. We use a least square minimization to determine a single gain factor which will best fit out model to our data. We plot the residuals to see if it has any significant structure. The null test is performed by subtracting corresponding data from one ref with the other. For Example, in the plot below, consider the first column. The first plot shows the NS slice of the rf0XX beam, with the NS slice of the rf1XX beam below it. The green points represent our data and errors, while the crimson curve represents the FEE model fitted to the data. The blue data points are the residuals between our data and the FEE model, while the orange curve is a 3rd order polynomial fit to it. The last plot displays the difference between the two plots above it. 

<p float="left">
  <img src="./imgs/null_test_XX_slices.png" width="49%" />
  <img src="./imgs/null_test_YY_slices.png" width="49%" />
</p>


#### Compare Beams

In this section we compare our tile maps with the FEE beam models of the MWA tiles using `compare_beams.py`. EW & NS slices of the tile maps, at the 4 pointings, are fit to corresponding slices of the FEE maps with a global gain factor obtained by least-squares minimization. A difference map is also plotted to show gradiants across the ratio of tile to FEE maps.

```
python compare_beams.py
```

<p float="left">
  <img src="./imgs/S35XX_rf1XX_0_beam_slices.png" width="49%" />
  <img src="./imgs/S35YY_rf1YY_0_beam_slices.png" width="49%" />
</p>


