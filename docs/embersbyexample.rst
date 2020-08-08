
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

