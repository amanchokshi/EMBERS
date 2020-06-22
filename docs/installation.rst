
=====================
 Installing Embers
=====================

.. currentmodule:: embers

Embers is designed to install cleanly with a single invocation
of the standard Python package tool

.. code-block:: console

    pip install embers

This should install embers and the small collection of Python libraries that embers depends on.
Using a virtual environment is suggested as it prevents errors from conflicting dependancies

.. code-block:: console

    python -m venv embers
    source emebrs/bin/activate
    pip install embers

Embers can also be installed manually from it's github repository

.. code-block:: console

    git clone https://github.com/amanchokshi/EMBERS.git
    cd EMBERS
    python -m venv embers
    source embers/bin/activate
    pip install .


If you find any problems or would like to suggest an improvement,
simply create an issue on the project’s GitHub page:

    https://github.com/amanchokshi/EMBERS

Good luck!
