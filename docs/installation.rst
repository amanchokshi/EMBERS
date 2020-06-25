
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

Conda virtual environment
--------------------------
EMBERS supports Python 3.6: 3.9. Using a conda environment is the easiest way to make sure that the correct version of Python is used. Begin by installing either `Anaconda <https://docs.anaconda.com/anaconda/install/>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. A conda environment called *embers* can be created as follows

.. code-block:: console

    conda create --name embers python=3.7
    conda activate embers
    pip install embers


Python virtual environment
--------------------------
Alternately, if you already have a correct version of Python installed, you can create a `venv <https://docs.python.org/3/library/venv.html/>`_ called *embers* using the following code, within which EMBERS and it's dependancies are cleanly installed

.. code-block:: console

    python -m venv embers
    source emebrs/bin/activate
    pip install embers



If you find any problems or would like to suggest an improvement,
simply create an issue on the project’s GitHub page:

    https://github.com/amanchokshi/EMBERS

Good luck!
