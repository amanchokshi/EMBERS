
===================
 Embers by Example
===================

The primary use case for *EMBERS* is the analysis of Radio Data obtained during the measurement of `MWA <https://www.mwatelescope.org/>`_ antenna beam-patterns
using satellites. The following examples outline the pipline used in the analysis of data for the Paper [Coming soon].

Setup
-----
*EMBERS* works best in a virtual environment and requires Python > 3.6. Follow the :doc:`Installation Instructions <installation>` to setup a suitable enviroment.

.. code-block:: console

    $ mkdir ~/EMBERS
    $ cd ~/EMBERS

The above creates an :samp:`EMBERS` directory, within which we can run our example code. Outputs will by default, be saved to :samp:`~/EMBERS/embers_out`.

A set of cli-tools enable easy interactions with the *EMBERS* packages, each of which can be executed with either of the following methods

.. code-block:: console

    $ python -m test_tool
    $ test_tool

These tools come built with help functions, which can be accessed at any time with the :samp:`--help` flag

.. code-block:: console

    $ test_tool --help

Raw Data Tree
^^^^^^^^^^^^^
It is important to know the directory structure of input data that *EMBERS* prefers. Data is organised into a root directory called :samp:`tiles_data` within which
sub directories for each MWA and reference tile exist (:samp:`S06XX`, :samp:`S06YY`, ......, :samp:`rf1XX`, :samp:`rf1YY`). Within each of these tile directories,
there exists a date directory for every day of observations in the :samp:`YYYY-MM-DD` format, within which live the raw RF data, in binary :samp:`.txt` files
saved every 30 minutes, with the naming convention as follows :samp:`S06XX_YYYY-DD-MM-hh:mm.txt`

.. code-block:: console

    tiles_data
    ├── S06XX
    │   ├── 2019-10-01
    │   │   ├── S06XX_2019-10-01-00:00.txt
    │   │   ├──         ........
    │   │   └── S06XX_2019-10-01-23:30.txt
    │   └── 2019-10-02
    │       ├── S06XX_2019-10-02-00:00.txt
    │       ├──         ........
    │       └── S06XX_2019-10-02-23:30.txt
    └── S06XX
        ├── 2019-10-01
        │   ├── S06XX_2019-10-01-00:00.txt
        │   ├──         ........
        │   └── S06XX_2019-10-01-23:30.txt
        └── 2019-10-02
            ├── S06XX_2019-10-02-00:00.txt
            ├──         ........
            └── S06XX_2019-10-02-23:30.txt

RF Tools
--------

:mod:`embers.rf_tools` is used to pre process, condition and preview raw rf data. Outputs of this module are saved to the :samp:`./embers_out/rf_tools` directory.

Waterfall Plots
^^^^^^^^^^^^^^^
To get a quick preview of the raw RF data, we create waterfall plots. The following code creates a waterfall plot of sample data provided with *EMBERS*

.. code-block:: console

    $ waterfall_single
    --------------------------------------------------
    No input data provided, using packaged sample data
    >>> waterfall_single --help, for more options
    --------------------------------------------------
    Waterfall plot saved to ./embers_out/rf_tools/S06XX_2019-10-10-02:30.png

.. image:: _static/imgs/waterfall_sample.png
    :width: 100%
    :alt: Waterfall Plot

We can also create a set of waterfall plots for all rf_files within a date interval

.. code-block:: console

    $ waterfall_batch --start_date=YYYY-MM-DD --stop_date=YYYY-MM-DD --data_dir=./tiles_data

Colormaps
^^^^^^^^^
*EMBERS* comes with two beautiful custom colormaps called :samp:`spectral` & :samp:`jade`. The :samp:`spectral` colormap is non-linear and is just used to 
visualise raw data and maximize dynamic range, while :samp:`jade` is perceptually uniform and sequential and is suitable for science. 
To get a preview of how amazing they are

.. code-block:: console

    $ colormaps

.. image:: _static/imgs/colormaps.png
    :width: 100%
    :alt: EMBERS custom colormaps

Align Data
^^^^^^^^^^

The RF Explorers used to record satellite data did not record data at exactly the same frequency and do not start recording at exactly the same time. 
In fact, the older models record at approximately 6 Hz, while the newer ones are capable of a sampling rate of nearly 9 Hz. This discrepency in sampling
rates makes it difficult to compare any two data samples. This issue is overcome by smoothing the data, along the time axis, with a Savitzky-Golay filter. 
Interpolating the smoothed data and resampling it at a constant frequency [ 0.5 Hz ] gives us a easier data set to work with. 

Two level of savgol filters are applied, first to capture deep nulls + small structure, and second level to smooth over noise. :samp:`align_single` can be
used to play with the various parameters available. Sensible defaults are provided as a starting point. The following code plots one frequency channel of
RF data and shows the efficacy of the selected smoothing filter.

.. code-block:: console
    
    $ align_single
    ----------------------------------------------------------
    No ref_file provided, using packaged sample data
    No tile_file provided, using packaged sample data
    No frequency channel provided, using 59 for sample data
    
    >>> savgol_interp_sample --help, for more options
    ----------------------------------------------------------
    Saving sample savgol_interp plot to: ./embers_out/rf_tools


.. image:: _static/imgs/align_data.png
    :width: 100%
    :alt: EMBERS custom colormaps

We can now align all the raw RF files within a date interval. Every pair of reference and MWA tile are smoothed and aligned and saved to compressed
:samp:`npz` file by :func:`~numpy.savez_compressed`.

.. code-block:: console

    $ align_batch --start_date=YYYY-MM-DD --stop_date=YYYY-MM-DD --data_dir=./tiles_data

