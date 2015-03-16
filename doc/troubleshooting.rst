.. _troubleshooting:

Troubleshooting
===============

Testing the Dependencies - Debugging
------------------------------------

Once installed, you should be able to import the python packages in a python console. 
MSNoise comes with a little script called `bugreport.py` that can be useful
to check if you have all the required packages (+ some extras).

The usage is such:

.. code-block:: sh

    $ msnoise bugreport -h

    usage: msnoise bugreport [-h] [-s] [-m] [-e] [-a]
    
    Helps determining what didn't work
    
    optional arguments:
      -h, --help     show this help message and exit
      -s, --sys      Outputs System info
      -m, --modules  Outputs Python Modules Presence/Version
      -e, --env      Outputs System Environment Variables
      -a, --all      Outputs all of the above


On my Windows machine, the execution of 

.. code-block:: sh

    $ msnoise bugreport -s -m

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


