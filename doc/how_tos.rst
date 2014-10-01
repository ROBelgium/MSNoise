.. warning:: This page is under construction, ...

.. _how_tos:

How To's
================

.. _multiproject:

Multi-Project
------------------------------------

Currently, MSNoise is NOT a regular Python Package, which means:

* It doesn’t have a setup.py file.
* It should not be installed in the lib/site-packages folder of the current python installation.
* It will create folder within its root folder (CROSS_CORRELATION, STACKS, MWCS, DTT folders). Thus:
* It should be installed in a folder that is writable to the user.

This should change in the future. Until then, when multiple users want to use MSNoise on the same machine, 
it would be better to have a copy of MSNoise for each (even: one copy for each "project").

Practically, for each "project", one has to download (or git fetch) the sources
of MSNoise in a new folder and run the installer there. If one's database
choice is MySQL, a new database must be created on the database server before
running the installer script. Depending on security recommandations from your
sysadmins, you might want to create specific users for specific databases. If
you are the only user and/or only running on your local machine, simply use the
same credentials for all databases.

The folder structure could look like this:

:: 

    .
    `-- MSNoiseProjects/
        |-- Project1/
        |   `-- STACKS/
        |   `-- s01scan_archive.py
        |   `-- ...
        `-- Project2/
        |   `-- STACKS/
        |   `-- s01scan_archive.py
        |   `-- ...
        `-- ..

Reprocess data
------------------------------------

When starting to use MSNoise, one will most probably need to re-run different
parts of the Workflow more than one time. By default, MSNoise is designed to
only process "what's new", which is antagonistic to what is wanted. Hereafter,
we present cases that will cover most of the re-run techniques:


When adding a new filter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If new filter are added to the filters list in the Configurator, one has to
reprocess all CC jobs, but not for filters already existing. The recipe is:

* Add a new filter, be sure to mark 'used'=1
* Set all other filters 'used' value to 0
* Redefine the flag of the CC jobs, from 'D'one to 'T'odo (see how)
* Run $python s03compute_cc.py
* Set back the other filters 'used' value to 1

The s03compute_cc.py will only compute the CC's for the new filter(s) and
output the results in the STACKS/ folder, in a subfolder named by a formatted
integer from the filter ID. For example: STACKS/01 for 'filter id'=1, STACKS/02
for 'filter id'=2, etc.


When changing the REF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO



When changing the MWCS parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the MWCS parameters are changed in the database, all MWCS jobs need to be
reprocessed.


TODO


When changing the dt/t parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TODO


Process workflow with more than 24 hours between steps
--------------------------------------------------------

TODO

Define one's own data structure of the waveform archive
---------------------------------------------------------

The data_structure.py file contains the known data archive formats.

.. note:: could this go to ObsPy somehow ?

TODO

How to have MSNoise work with 2+ data structures at the same time
-------------------------------------------------------------------

TODO

How to duplicate/dump the MSNoise configuration
------------------------------------------------

TODO

Get the plots from the MSNoise SRL paper
------------------------------------------

TODO



.. _testing:
Testing the Dependencies
------------------------------------

Once installed, you should be able to import the python packages in a python console. 
MSNoise comes with a little script called `bugreport.py` that can be useful
to check if you have all the required packages (+ some extras).

The usage is such:

.. code-block:: sh

    $ python bugreport.py -h

    usage: bugreport.py [-h] [-s] [-m] [-e] [-a]
    
    Helps determining what didn't work
    
    optional arguments:
      -h, --help     show this help message and exit
      -s, --sys      Outputs System info
      -m, --modules  Outputs Python Modules Presence/Version
      -e, --env      Outputs System Environment Variables
      -a, --all      Outputs all of the above


On my Windows machine, the execution of 

.. code-block:: sh

    $ python bugreport.py -s -m

results in:

.. code-block:: sh

    ************* Computer Report *************
    
    ----------------+SYSTEM+-------------------
    Windows
    seis31
    7
    6.1.7601
    AMD64
    Intel64 Family 6 Model 42 Stepping 7, GenuineIntel
    
    ----------------+PYTHON+-------------------
    Python: 2.7.5 |Anaconda 1.7.0 (64-bit)| (default, Jul  1 2013, 12:37:52) [MSC v.1500 64 bit (AMD64)]
    
    ---------------+MODULES+-------------------
    
    Required:
    [X] numpy: 1.7.1
    [X] scipy: 0.12.0
    [X] pandas: 0.12.0
    [X] matplotlib: 1.3.0
    [X] statsmodels: 0.5.0
    [X] sqlalchemy: 0.8.2
    [X] traitsui: 4.3.0
    [X] traits: 4.3.0
    [X] enable: 4.3.0
    [X] scikits.samplerate: present (no version)
    [X] obspy: present (no version)
    [X] sphinx: 1.1.3
    [X] jinja2: 2.7.1
    
    Backends: (at least one is required)
    [X] wx: 2.8.12.1
    [ ] PyQt4: not found
    [X] PySide: 1.2.1
    
    Not required, just checking:
    [X] setuptools: 0.6
    [X] reportlab:  $Id$
    [X] configobj: 4.7.2
    [X] pkg_resources: present (no version)
    [ ] paramiko: not found
    [X] ctypes: 1.1.0
    [X] pyparsing: 1.5.6
    [X] distutils: 2.7.5
    [X] IPython: 1.0.0
    [X] vtk: present (no version)

The [X] marks the presence of the module. In the case above, PyQt4 is missing, but that's not a problem because
`wx` or `PySide` are present, so traitsui has a backend to render the GUI for the Configurator. The "not-required"
packages are checked for information, those packages can be useful for reporting / hacking / rendering the data.

To install a missing package, for example *obspy*, use the easy_install command (easy_install is a python script that
comes with setuptools):

.. code-block:: sh

    $ easy_install obspy




