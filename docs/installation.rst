
=====================
 Installing Embers
=====================

.. currentmodule:: embers

Embers is designed to install cleanly with a single invocation
of the standard Python package tool

.. code-block:: console

    pip install embers

Embers can also be installed manually from it's github repository

.. code-block:: console

    git clone https://github.com/amanchokshi/EMBERS.git
    cd EMBERS
    pip install .

These methods should install embers and the small collection of Python libraries that embers depends on.
Using a virtual environment is suggested as it prevents errors from conflicting dependancies

The embers package can now be imported with :samp:`import embers`

Pyenv virtual environment
-------------------------
Ember supports Python 3.6-3.9. Using a :samp:`pyenv` environment is the cleanest route to the correct pthon version. Install pyenv with the instructions at `<https://realpython.com/intro-to-pyenv/>`_. Now, to install the correct version of python and make a virtual environment

.. code-block:: console
    
    # Install python version 3.8.3
    pyenv install 3.8.3

    # Create a virtual environment called embers
    # env located at ~/.pyenv/versions/embers
    pyenv virtualenv 3.8.3 embers

    # virtual env automatically activated when in EMBERS
    # Link the EMBERS dir with the virtual env
    # creates a .python-version file
    python local embers

    # Install the embers package & dependancies
    pip install embers


Conda virtual environment
--------------------------
EMBERS supports Python 3.6-3.9. Using a :samp:`conda` environment is the easiest way to make sure that the correct version of Python is used. Begin by installing either `Anaconda <https://docs.anaconda.com/anaconda/install/>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. A conda environment called *embers-env* can be created as follows

.. code-block:: console

    conda create --name embers-env python=3.7
    conda activate embers-env
    pip install embers


Python virtual environment
--------------------------
Alternately, if you already have a correct version of Python installed, you can create a `venv <https://docs.python.org/3/library/venv.html/>`_ called *embers-env* using the following code, within which EMBERS and it's dependancies are cleanly installed

.. code-block:: console

    python -m venv embers-env
    source embers-env/bin/activate
    pip install embers



If you find any problems or would like to suggest an improvement,
simply create an issue on the project’s GitHub page:

    https://github.com/amanchokshi/EMBERS

Good luck!
