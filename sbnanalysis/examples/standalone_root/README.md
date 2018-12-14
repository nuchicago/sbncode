Standalone ROOT Scripts
=======================
Author: A. Mastbaum <mastbaum@uchicago.edu>, 2018/12/14

This is an example of working with sbnanalysis output files using a standalone
compiled ROOT program.

Notes
-----
* The Makefile links to the shared library containing the ROOT dictionary for
  the Event class, `libsbnanalysis_Event.so`
* The script loops through the event tree and prints some event information

Building
--------
To build, set up ROOT and run:

    $ make

Running
-------
To run the example, build then run:

    $ ./eventdump some_sbnanalysis_output.root

