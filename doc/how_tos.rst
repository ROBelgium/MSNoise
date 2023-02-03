.. _how_tos:

How To's - Recipes
==================

Run MSNoise only to have Power Spectral Densities and Spectrograms:
-------------------------------------------------------------------

This recipe is a kind of "let's check noise content rapidly":

.. code-block:: sh

    msnoise db init --tech 1

    msnoise config set startdate=2019-01-01
    msnoise config set enddate=2019-02-01
    msnoise config set response_path=/path/to/response_files

    msnoise scan_archive --path /path/to/archive --recursively
    msnoise populate --fromDA
    msnoise db update_loc_chan
    msnoise new_jobs --init --nocc

    msnoise admin # check the qc_* parameters

    msnoise qc compute_psd
    msnoise qc plot_psd YA.UV05.00.HHZ


Run the simplest MSNoise run ever
---------------------------------

This recipe is a kind of "let's check this data rapidly":

.. code-block:: sh

    msnoise db init --tech 1

    msnoise config set startdate=2019-01-01
    msnoise config set enddate=2019-02-01
    msnoise config set overlap=0.5
    msnoise config set mov_stack=1,5,10

    msnoise scan_archive --path /path/to/archive --recursively
    msnoise populate --fromDA
    msnoise db update_loc_chan
    msnoise new_jobs --init

    msnoise admin # add 1 filter in the Filter table
    # or
    msnoise db execute "insert into filters (ref, low, mwcs_low, high, mwcs_high, mwcs_wlen, mwcs_step, used) values (1, 0.1, 0.1, 1.0, 1.0, 12.0, 4.0, 1)"

    msnoise cc compute_cc
    msnoise cc stack -r
    msnoise reset STACK
    msnoise cc stack -m
    msnoise cc dvv compute_mwcs
    msnoise cc dvv compute_dtt
    msnoise cc dvv plot dvv


Run MSNoise using lots of cores on a HPC
----------------------------------------

Avoid Database I/O by using the ``hpc`` flag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With MSNoise 1.6, most of the API calls have been cleaned from calling the
database, for example the ``def stack()`` called a SELECT on the database for
each call, which is useless as configuration parameters are not supposed to
change during the execution of the code. This modification allows running
MSNoise on an HPC infrastructure with a remote central MySQL database.

The new configuration parameter ``hpc`` is used for flagging if MSNoise is 
running High Performance. If True, the jobs processed at each step are marked
Done when finished, but the next jobtype according to the workflow is not
created. This removes a lot of select/update/insert actions on the database
and makes the whole much faster (one INSERT instead of tons of 
SELECT/UPDATE/INSERT).

Commands and actions with ``hpc`` = N :

* ``msnoise new_jobs``: creates the CC jobs
* ``msnoise cc compute_cc``: processes the CC jobs and creates the STACK jobs
* ``msnoise cc stack -m``: processes the STACK jobs and creates the MWCS jobs
* etc...

Commands and actions with ``hpc`` = Y :

* ``msnoise new_jobs``: creates the CC jobs
* ``msnoise cc compute_cc``: processes the CC jobs
* ``msnoise new_jobs --hpc CC:STACK``: creates the STACK jobs based on the CC 
  jobs marked "D"one
* ``msnoise cc stack -m``: processes the STACK jobs
* ``msnoise new_jobs --hpc STACK:MWCS``: creates the MWCS jobs based on the 
  STACK jobs marked "D"one
* etc...

Set up the HPC
~~~~~~~~~~~~~~

To avoid having to rewrite MSNoise for using techniques relying on MPI or other
parallel computing tools, I decided to go "simple", and this actually works. The
only limitation of the following is that you need to have a strong MySQL server
machine that accepts hundreds or thousands of connections. In my case, the
MySQL server is running on a computing blade, and its my.cnf is configured to
allow 1000 users/connections, and to listen on all its IPs.

The easiest set up (maybe not your sysadmin's preferred, please check), is to

* install miniconda on your home directory and make miniconda's python
  executable your default python (I add the paths to .profile).
* Then install the requirements and finally MSNoise.
* As usual, create a project folder and ``msnoise db init`` there, choose MySQL
  and provide the hostname of the machine running the MySQL server.

At that point, your project is ready. I usually request an interactive node on
the HPC for doing the ``msnoise populate`` and ```msnoise scan_archive``. Our
jobs scheduler is PBS, so this command

.. code-block:: sh

    qsub -I -l walltime=02:00:00 -l select=1:ncpus=16:mem=1g

requests an Interactive node with 16 cpus, 1GB ram, for 2 hours. Once connected,
check that the python version is correct (or source .profile again). Because
we requested 16 cores, we can ``msnoise -t 16 scan_archive --init``.

Depending on the server configuration, you can maybe run the ``msnoise admin``
on the login node, and access it via its hostname:5000 in your browser. If not,
the easiest way to set up the config is running
``msnoise config set <parameter>=<value>`` from the console. To add filters,
do it either:

* in the Admin 
* using MySQL workbench connected to your MySQL server
* using such commands ``msnoise db execute "insert into filters (ref, low, mwcs_low, high, mwcs_high, mwcs_wlen, mwcs_step, used) values (1, 0.1, 0.1, 1.0, 1.0, 12.0, 4.0, 1)"``
* using ``msnoise db dump``, edit the filter table in CSV format, then ``msnoise db import filters --force``

Once done, the project is set up and should run. Again, test if all goes OK in
an interactive node.

To run on N cores in parallel, we have the advantage that, e.g. for CC jobs, the
day-jobs are independent. We can thus request an "Array" of single cores, which
is usually quite easy to get on HPCs (most users run heavily parallel codes and
request large number of "connected" cores, while we can run "shared").

The job file in my PBS case looks like this for computing the CC:

.. code-block:: sh

    #!/bin/bash
    #PBS -N MSNoise_PDF_CC
    #PBS -l walltime=01:00:00
    #PBS -l select=1:ncpus=1:mem=1g
    #PBS -l place=shared
    #PBS -J 1-400
    cd /scratch-a/thomas/2019_PDF
    source /space/hpc-home/thomas/.profile
    msnoise cc compute_cc

This requests 400 cores with 1GB of RAM. The content of my .profile file
contains:

.. code-block:: text

    # added by Miniconda3 installer
    export PATH="/home/thomas/miniconda3/bin:$PATH"
    export MPLBACKEND="Agg"

The last line is important as nodes are usually "head-less" and matplotlib and
packages relating to it would fail if they expect a gui-capable system.

For submitting this job, run ``qsub qc.job``. The process usually routes stdout
and stderr to files in the current directory, make sure to check them if jobs
seem to have failed. If all goes well, calling ``msnoise info -j`` repeatedly
from the login or interactive node's console should show the evolution of Todo,
In Progress and Done jobs.

.. note:: HPC experts are welcome to suggest, comment, etc... It's a quick'n'dirty
    solution, but it works for me!


Reprocess data
--------------

When starting to use MSNoise, one will most probably need to re-run different
parts of the Workflow more than one time. By default, MSNoise is designed to
only process "what's new", which is antagonistic to what is wanted. Hereafter,
we present cases that will cover most of the re-run techniques:


When adding a new filter
~~~~~~~~~~~~~~~~~~~~~~~~

If new filter are added to the filters list in the Configurator, one has to
reprocess all CC jobs, but not for filters already existing. The recipe is:

* Add a new filter, be sure to mark 'used'=1
* Set all other filters 'used' value to 0
* Redefine the flag of the CC jobs, from 'D'one to 'T'odo with the following:
* Run ``msnoise reset CC --all``
* Run ``msnoise cc compute_cc``
* Run next commands if needed (stack, mwcs, dtt)
* Set back the other filters 'used' value to 1

The compute_cc will only compute the CC's for the new filter(s) and
output the results in the STACKS/ folder, in a sub-folder named by a formatted
integer from the filter ID. For example: STACKS/01 for 'filter id'=1, STACKS/02
for 'filter id'=2, etc.


When changing the REF
~~~~~~~~~~~~~~~~~~~~~

When changing the REF (``ref_begin`` and ``ref_end``), the REF stack has to be
re-computed:

.. code-block:: sh

    msnoise reset STACK --all
    msnoise cc stack -r

The REF will then be re-output, and you probably should reset the MWCS jobs to
recompute daily correlations against this new ref:

.. code-block:: sh

    msnoise reset MWCS --all
    msnoise cc dvv compute_mwcs


When changing the MWCS parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the MWCS parameters are changed in the database, all MWCS jobs need to be
reprocessed:

.. code-block:: sh

    msnoise reset MWCS --all
    msnoise cc dvv compute_mwcs

shoud do the trick.


When changing the dt/t parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: sh

    msnoise reset DTT --all
    msnoise cc dvv compute_dtt


Recompute only the specific days
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You want to recompute CC jobs after a certain date only, for whatever reason:

.. code-block:: sh

    msnoise reset CC --rule="day>='2019-01-01'"

SQL experts can also use the ``msnoise db execute`` command (with caution!):

.. code-block:: sh

    msnoise db execute "update jobs set flag='T' where jobtype='CC' and day>='2019-01-01'"

If you want to only reprocess one day:

.. code-block:: sh

    msnoise reset CC --rule="day='2019-01-15'"



Define one's own data structure of the waveform archive
-------------------------------------------------------

The data_structure.py file contains the known data archive formats. If another
data format needs to be defined, it will be done in the ``custom.py`` file
in the current project folder:

.. seealso:: Check the "Populate Station Table" step in the :doc:`workflow/002_populate`.


How to have MSNoise work with 2+ data structures at the same time
-----------------------------------------------------------------

In this case, the easiest solution is to scan the archive(s) with the "Lazy
Mode":

.. code-block:: sh

    msnoise scan_archive --path /path/to/archive1/ --recursively
    msnoise scan_archive --path /path/to/archive2/ --recursively

etc.

Remember to either manually fill in the station table, or

.. code-block:: sh

    msnoise populate --fromDA



How to duplicate/dump the MSNoise configuration
-----------------------------------------------

To export all tables of the current database, run

.. code-block:: sh

    msnoise db dump

This will create as many CSV files as there are tables in the database.

Then, on a new location, init a new msnoise project and import the tables
one by one:

.. code-block:: sh

    msnoise db init
    msnoise db import config --force
    msnoise db import stations --force
    msnoise db import filters --force
    msnoise db import data_availability --force
    msnoise db import jobs --force


Check if my response file works
-------------------------------

To check if your response file can be used by msnoise, you simply should
check that it is readable with ObsPy and contains the response information.
In a python shell, do the following:

.. code-block:: python

    from obspy.core import UTCDateTime, read_inventory, read
    st = read("/path/to/a/file/for/station/XX.BBB")
    inv = read_inventory("/path/to/the/response/for/station/XX.BBB)
    print(inv)
    response = inv.get_response(st[0].id, st[0].stats.starttime)
    print(response)


alternatively, if you have configured the path to the response files
(``response_path``) correctly you can also call the msnoise api:

.. code-block:: python

    from msnoise.api import connect, preload_instrument_responses
    st = read("/path/to/a/file/for/station/XX.BBB")
    db = connect()
    inv = preload_instrument_responses(db, return_format="inventory")
    response = inv.get_response(st[0].id, st[0].stats.starttime)
    print(response)

.. _testing:

Testing the Dependencies
------------------------

Once installed, you should be able to import the python packages in a python console. 
MSNoise comes with a little script called `bugreport.py` that can be useful
to check if you have all the required packages (+ some extras).

The usage is such:

.. code-block:: sh

    $ msnoise bugreport -h

    usage: msnoise bugreport [-h] [-s] [-m] [-e] [-a]
    
    Helps determining what didn\'t work
    
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
    PC1577-as
    10
    10.0.17134
    AMD64
    Intel64 Family 6 Model 158 Stepping 9, GenuineIntel
    
    ----------------+PYTHON+-------------------
    Python:3.7.3 | packaged by conda-forge | (default, Jul  1 2019, 22:01:29) [MSC v.1900 64 bit (AMD64)]
    
    This script is at d:\pythonforsource\msnoise_stack\msnoise\msnoise\bugreport.py
    
    ---------------+MODULES+-------------------
    
    Required:
    [X] setuptools: 41.2.0
    [X] numpy: 1.15.4
    [X] scipy: 1.3.0
    [X] pandas: 0.25.0
    [X] matplotlib: 3.1.1
    [X] sqlalchemy: 1.3.8
    [X] obspy: 1.1.0
    [X] click: 7.0
    [X] pymysql: 0.9.3
    [X] flask: 1.1.1
    [X] flask_admin: 1.5.3
    [X] markdown: 3.1.1
    [X] wtforms: 2.2.1
    [X] folium: 0.10.0
    [X] jinja2: 2.10.1
    
    Only necessary if you plan to build the doc locally:
    [X] sphinx: 2.2.0
    [X] sphinx_bootstrap_theme: 0.7.1
    
    Graphical Backends: (at least one is required)
    [ ] wx: not found
    [ ] pyqt: not found
    [ ] PyQt4: not found
    [X] PyQt5: present (no version)
    [ ] PySide: not found
    
    Not required, just checking:
    [X] json: 2.0.9
    [X] psutil: 5.6.3
    [ ] reportlab: not found
    [ ] configobj: not found
    [X] pkg_resources: present (no version)
    [ ] paramiko: not found
    [X] ctypes: 1.1.0
    [X] pyparsing: 2.4.2
    [X] distutils: 3.7.3
    [X] IPython: 7.7.0
    [ ] vtk: not found
    [ ] enable: not found
    [ ] traitsui: not found
    [ ] traits: not found
    [ ] scikits.samplerate: not found


The [X] marks the presence of the module. In the case above, PyQt4 is missing, but that's not a problem because
`PyQt5` is present. The "not-required" packages are checked for information, those packages can be useful for reporting / hacking / rendering the data.


