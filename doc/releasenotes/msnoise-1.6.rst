.. include:: ../configs.hrst

MSNoise 1.6
===========

Release date: XX August 2019


Release type: major


Release notes:

.. contents::
    :local:

Introduction
------------
More than 2 years after the last major release (:doc:`msnoise-1.5`) I'm proud
to announce the new :doc:`msnoise-1.6`. It is a **major** release, with a
massive amount of work since the last one: in `GitHub numbers
<https://github.com/ROBelgium/MSNoise/graphs/contributors?from=2017-04-28&to
=2019-08-12&type=c>`_
, it's over TODO commits and over TODO lines of code and documentation changed
or added!

*End of summer is also a very special period for MSNoise, as it has been 9
years since Corentin contacted Florent and that I immediately started working
on this package. 2010-2019. Nine years. Wow. MSNoise has now a few thousand
lines of code and more than 100 pages of documentation, it is widely used and
scientists around the globe use it and even make super cool publications out 
of their results! So proud!* 


MSNoise 1.6 introduces a series of **new features** :

* The workflow has been rewritten to create "job types" per step, making it 
  easier for users to reset a jobs before a specific step. 
* This and other smaller adaptations to the code allows to run MSNoise more
  effiently, e.g. on a HPC (see hpc_).
* The components to compute can be defined for "single-station" and 
  "cross-station" independently.
* The compute_cc step has been completely rewritten to be much, much faster.
* A new ``db`` top level command, which among others include ``dump`` and
  ``import`` commands for handling the tables from the database.


.. todo

As always, this version has benefited from outputs/ideas/pull 
requests/questions from several users/friends.

Thanks to all for using MSNoise, and please, let us know why/how you use it
(and please cite it!)!

To date, we found/are aware of 70 publications using MSNoise! That's the best
validation of our project ever and it has doubled since last release!! 


*Thomas*


~~~~

PS: if you use MSNoise for your research and prepare publications, **please
consider citing it**:

**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.


Requirements
------------

This version will be the last to be tested on Python 2.7. The EOL
(end of life) of 2.7 is 2020, which means it is high time for users to migrate.
For users having a complete set of tools in Python 2.7 and not keen to move 
to 3.x soon, the incredible easiness of creating a Python 3.x environment in 
conda, for example, will allow them to run MSNoise in the future.

There were no changes in the requirements. Note that MSNoise is always tested
against the latest release versions of the main packages, so older installations
that are not maintained/updated regularly (years) could encounter issues. 
Please make sure you have the latest version of Numpy and Scipy (and MKL), as
performance gets better and better (especially since Anaconda Inc. released 
its fast MKL implementations for all users, in the conda-forge channel). 



Configuration Parameters
------------------------

* ADDED: ``hpc`` for flagging if MSNoise is running High Performance. If 
  True, the jobs processed at each step are marked Done when finished, but the
  next jobtype according to the workflow is not created. This removes a lot 
  of select/update/insert actions on the database and makes the whole much 
  faster (See hpc_).

* ADDED: ``archive_format`` will tell the `obspy.core.stream.read ObsPy
  function
  <https://docs.obspy.org/packages/autogen/obspy.core.stream.read.html>`_ the
  format of the files to read in the archive during the `scan_archive` stage.
  If left empty, `obspy.core.stream.read` will automatically detect the format
  of the file, which results in a slightly slower reading.

* ADDED: ``whitening_type`` to allow for the standard brutal whitening (all-to-
  1.0) or whitening by dividing by the PSD of the singnal.

* ADDED: ``components_to_compute_single_station`` to separate the components
  computed for single stations and cross-stations. For example, the ``ZZ`` can
  be computed for STA1 vs STA2 and, at the same time, ``EZ,EN,EZ`` for STA1 and
  for STA2.

* CHANGED: ``cronday`` should now be a positive float number (negative numbers
  are still accepted for backward compatibility) or a string designating any
  optional number of weeks, days and/or hours in the format 'Xw Xd Xw' in this
  order, where a week represents 7 days (some valid examples would be:
  '1w 3d 6h', '2w 1d', '1d 12h', or '6h').  Its meaning is unchanged: if the
  ``--init`` option is not used, `scan_archive` will only scan files modified
  during the time in the past designated by `crondays` .


Workflow changes
----------------
The workflow has been rewritten to create "job types" per step, making it 
easier for users to reset a jobs before a specific step.

New job types are:
  
* A ``STACK`` job is created when a ``CC`` job is successful
* A ``MWCS`` job is created when a ``STACK`` job is successful
* A ``DTT`` job is created when a ``MWCS`` job is successful

About the ``STACK`` jobs, it is important to call the "ref stack" before the
"mov stack", as the "ref stack" will run on the ``STACK`` jobs, check if 
any date matches the date range of the ``ref_begin`` - ``ref_end``, do the 
stacks if needed, then will **not** reset the ``STACK`` jobs to ``Todo`` so 
that "mov stacks" can be done. 

TODO:: Note that calling ``msnoise stack -r -m`` will 
execute the steps in the right order.



Preprocessing and Cross-Correlation
-----------------------------------

Preprocessing
~~~~~~~~~~~~~

* Only small changes were done for this step, mainly checks of matching
  sampling rates, empty streams. DB-related optimisations make this step 
  faster too.

Cross-Correlation TODO
~~~~~~~~~~~~~~~~~~~~~~

* The ``compute_cc`` step has been completely rewritten to make use of 2D arrays
  holding the data, processing them "in place" for the different steps (FFT,
  whitening, etc). This results in much more efficient computation. The process
  slides on time windows and computes the correlations using indexes in a 2D
  array, therefore avoiding an exponential number of identical operations on
  data windows.

* This new code is the default ``compute_cc``, and it doesn't allow computing
  rotated components. For users needing ``R`` or ``T`` components, there are two
  options: either use the old code, now named ``compute_cc_rot``, or compute the
  full (6 components actually are enough) tensor using the new code, and rotate
  the components afterwards. From initial tests, this latter solution is a lot
  faster than the first, thanks to the new processing in 2D.

* It is now possible to do the Cross-Correlation (classic "CC"), the Auto-
  Correlation ("AC") or the Cross-Components within the same station ("SC").
  To achieve this, we removed the `ZZ`, `ZT`, etc parameters from the
  configuration and replaced it with ``components_to_compute`` which takes a list:
  e.g. `ZZ,ZE,ZN,EZ,EE,EN,NZ,NE,NN` for the full non-rotated tensor between
  stations. Adding components to the new
  ``components_to_compute_single_station`` will allow computing the
  cross-components (SC) or auto-correlation (AC) of each station.

* The cross-correlation is done on sliding windows on the available data. For
  each window, if one trace contains a gap, it is eliminated from the
  computation. This corrects previous errors linked with gaps synchronised in
  time that lead to perfect sinc autocorrelation functions. The windows should
  have a duration of at least "2 times the `maxlag`+1" to be computable.

.. todo:: params.whiten is not used in the new compute_cc

* The whitening procedure can be skipped by setting the ``whitening``
  configuration to `None`. The two other ``whitening`` modes are "[A]ll except
  for auto-correlation" or "Only if [C]omponents are different". This allows
  skipping the whitening when, for example, computing ZZ components for very
  close by stations (much closer than the wavelength sampled), leading to
  spatial autocorrelation issues. "A" and "N" are supported by the new
  ``compute_cc``, while "C" is only supported by the old version, now named
  ``compute_cc_rot``.


Command Line changes
--------------------


Top level DB command
~~~~~~~~~~~~~~~~~~~~

I've added a new command group called ``db`` that gathers all db-related actions:

* ``msnoise db init`` is a replacement for the ``msnoise install``
* ``msnoise db upgrade`` is a replacement for the ``msnoise upgrade_db``
* ``msnoise db clean_duplicates`` deletes duplicate jobs (might happen). Unique
  sets of ``day``, ``pair`` and ``jobtypes`` are considered.
* ``msnoise db execute`` allows executing SQL queries on the database (Expert 
  Mode). If the query includes a SELECT command, the result is returned as a
  table.
* ``msnoise db dump`` allows dumping all tables from the database to CSV files
* ``msnoise db import`` allows importing individual tables and replace them in
  the database.


The ``config`` command group has been reworked and the ``get`` subcommand has
been added to retrieve the values of a list of configuration parameters:

* ``msnoise config get <param>`` will display the value of the configuration
  parameter `<param>`.
* ``msnoise config set <param>=<value>`` will set the value of the configuration
  parameter `<param>` to `<value>`.
* ``msnoise config gui`` will run the deprecated configuration graphical
  interface (that was previously available through ``msnoise config`` but
  supersed by the web interface using ``msnoise admin``).
* ``msnoise config sync`` synchronises station metadata from inventory/dataless
  (was previously available as ``msnoise config --sync``).


Other changes
~~~~~~~~~~~~~

* ``msnoise info`` also prints information stored in the `db.ini` file (for
  security reason, the password is masked in the output but be aware that it is
  still stored in clear text in the file.).
* ``msnoise info -j`` reports all jobs types, including those of plugins.
* Added the possibility to walk in subfolders recursively by using 
  ``scan_archive --path -r``

Note, all commands are documented: :doc:`../clickhelp/msnoise`.


API Changes
-----------

* New ``get_params`` API method to return a ``Params`` class containing all 
  configuration bits, avoiding unnecessary calls to the DB later. This method
  will automatically populate from the ``defaults`` and return elements 
  having the right type. 
* Many (many) small optimizations to the core functions to make them faster,
  mostly when interacting with the DB, thanks to ``get_params``.
* ``update_config`` can now modify parameters for installed plugins.
* New ``massive_update_job`` to update a list of ``Jobs`` to a given flag.
* ``build_ref_datelist`` and ``build_movstack_datelist`` return the smallest 
  date from the data_availability table if the ``ref_begin`` or ``startdate``
  have not been modified from their default value of ``1970-01-01``.
* Removed ``linear_regression`` as this is now included in ObsPy.
* Modified ``get_dtt_next_job`` returns jobs in random order.
* Added missing documentation for several methods.

See :doc:`../api`.


Performance and Code improvements
---------------------------------

.. _hpc:

High Performance - Reducing DB access
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most of the API calls have been cleaned from calling the database, for 
example the ``def stack()`` called a SELECT on the database for each call, 
which is useless as configuration parameters are not supposed to change 
during the execution of the code. This modification allows running MSNoise on
an HPC infrastructure with a remote central MySQL database.

The new configuration parameter ``hpc`` is used for flagging if MSNoise is 
running High Performance. If True, the jobs processed at each step are marked
Done when finished, but the next jobtype according to the workflow is not
created. This removes a lot of select/update/insert actions on the database
and makes the whole much faster (one INSERT instead of tons of 
SELECT/UPDATE/INSERT).

Commands and actions with ``hpc`` = N :

* ``msnoise new_jobs``: creates the CC jobs
* ``msnoise compute_cc``: processes the CC jobs and creates the STACK jobs
* ``msnoise stack -m``: processes the STACK jobs and creates the MWCS jobs

Commands and actions with ``hpc`` = Y :

* ``msnoise new_jobs``: creates the CC jobs
* ``msnoise compute_cc``: processes the CC jobs
* ``msnoise new_jobs --hpc CC:STACK``: creates the STACK jobs based on the CC 
  jobs marked "D"one
* ``msnoise stack -m``: processes the STACK jobs
* ``msnoise new_jobs --hpc STACK:MWCS``: creates the MWCS jobs based on the 
  STACK jobs marked "D"one

.. _scan-archive:

Rework of scan_archive
~~~~~~~~~~~~~~~~~~~~~~

The code behind the ``scan_archive`` command has been deeply reworked to ease
maintenance and debugging, and the reading of the files and directories has
been improved thanks to the `scandir function
<https://github.com/benhoyt/scandir>`_ included in Python 3.5 and backported in
the `scandir` module.  The multiprocessing strategy (used if the
``-t``/``--threads`` option is provided) has also been reworked to limit the
expensive step of process creation: the list of directories to scan is now
split at startup and each child process concurrently treats an equal part of
it.

.. _table-prefix:

Prefix of database tables
~~~~~~~~~~~~~~~~~~~~~~~~~

The new ``msnoise db init`` project initialisation command now prompts you to
choose an optional prefix for the name of the tables that will be created in
the chosen database.  If you enter the prefix `myprefix`, the tables will be
named using the prefix `myprefix_`.  This allows to share the same database
with several MSNoise projects.

Comparison with MSNoise 1.2 as published in the SRL article:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
In the 2014 SRL article, we computed the dv/v for the UnderVolc project on a 
4vCPU virtual machine running on a powerful ESX system. For this test, I 
analyzed the same data set on a 4 year-old 16 CPU blade. The timings 
mentioned below are then multiplied by 4 to account for the CPU number 
difference.

Summing the total time needed, we reach 37 hours for the SRL version, and 12 
hours for MSNoise 1.6. The seedup is not fully linear, as the current code 
supports running on over 500 CPUs (as long as your MySQL server can handle 
it) but no MySQL server could have handled the so many connections/requests 
from the old version. The ``compute_cc2`` computation time scales roughly 
linearly with the amount of components, contrary to the old ``compute_cc`` 
which was exponential).

	
+--------------+--------------+---------------+---------------+--------+
|STEP          |v1.6 (16 CPU) |v1.6 (4 CPU)   |SRL (4 vCPU)   |SPEEDUP |
+==============+==============+===============+===============+========+
|scan_archive  |385 seconds   |1540 seconds   |1800 seconds   |1.2x    |
+--------------+--------------+---------------+---------------+--------+
|new_jobs      |27 seconds    |27 seconds     |1800 seconds   |66.0x   |
+--------------+--------------+---------------+---------------+--------+
|compute_cc2   |4817 seconds  |19268 seconds  |75600 seconds  |3.9x    |
+--------------+--------------+---------------+---------------+--------+
|stack -r      |58 seconds    |232 seconds    |1980 seconds   |8.5x    |
+--------------+--------------+---------------+---------------+--------+
|stack -m      |1124 seconds  |4496 seconds   |21600 seconds  |4.8x    |
+--------------+--------------+---------------+---------------+--------+
|compute_mwcs  |4209 seconds  |16836 seconds  |28800 seconds  |1.7x    |
+--------------+--------------+---------------+---------------+--------+
|compute_dtt   |264 seconds   |1056 seconds   |3600 seconds   |3.4x    |
+--------------+--------------+---------------+---------------+--------+
|Total         |10884 seconds |43455 seconds  |135180 seconds |3.1x    |
+--------------+--------------+---------------+---------------+--------+
|Total (hours) |3 hours       |12 hours       |37 hours       |3.1x    |
+--------------+--------------+---------------+---------------+--------+

Next steps of improving this workflow will be:

* Reduce I/O: by storing the ``mov_stack`` elements in two-dimensional HDF5 
  dataframes
* Reduce CPU: by computing the ``mov_stack`` in a 2D Pandas DataFrame directly
* Reduce CPU: by pre-computing MWCS windows and pair-wise computing the 
  delays (reduces drastically the number of FFT calls)

**If anyone feels like focusing on those aspects and providing Pull Requests, 
welcome!**


Other changes
-------------

Web-based Admin Interface Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Feature: The pagination size (100, 200, etc. row) is now allowed on the
  Station, Job and DataAvailability views.
* Feature: The Config list view can be sorted by name.
* Feature: The home page has changed to show all job types and exposes 
  buttons to execute actions equivalent to `msnoise reset` in the console.

Other Bugfixes TODO
~~~~~~~~~~~~~~~~~~~

* Removed the call to ``scipy.stats.nanmean`` and replaced by ``numpy.nanmean``
* Better error message in compute_cc when the content of the slice is only zeros
  or smaller than ``rms_threshold``
* Checked all "integer" - "float" warnings from numpy/scipy
* crondays were hardcoded to -2, now taking the ``crondays`` value from the DB
* Py3 error when ``msnoise scan_archive`` in cron mode
* ``msnoise info`` does not crash any more if some station information are
  `NULL` in the database


Plot Updates TODO
~~~~~~~~~~~~~~~~~
Not sure we did modify something here

See :doc:`../plotting`.

Documentation TODO
~~~~~~~~~~~~~~~~~~

* New elements for configuring MySQL and MariaDB, thanks to Lukas Preiswerk.
* The description of the new steps and the HPC mode.



Upgrading an existing project to MSNoise 1.6
--------------------------------------------

Some users will want to keep their current project without recomputing
everything. This requires adding a few configuration parameters to the database

Running the following command will take care of the upgrade from 1.5 to 1.6:

.. code-block:: sh

    msnoise db upgrade


A final note about development pace and choices
-----------------------------------------------

* MSNoise team is

  * **1 developper** (Thomas)
  * 1 dedicated debugger (Corentin)
  * less than 10 really *active* users, providing feedback and/or lines of codes
    (Esteban, Raphaël, Aurélien, Carmelo, Clare, Rob ...)

* All software engineering ideas are coming from too infrequent beerstormings
  between Thomas & others
* The web-interface and the plugin support were developed during Thomas'
  holidays

If you need help, please ask your questions on the mailing list. Don't be afraid
to ask. If you have ideas, please share them. If you develop codes to supplement
MSNoise, please share them, even if very small, even if you don't master gitHub.
If you have complaints, post them too, but remember that the package you are
using has been coded by 1 person, and that it's not his full time job. So
MSNoise is provided "as-is", carefully written and tested, but there will be
bugs, issues, incompatibility with certain python installations, OS or module
versions. If you **want or need** developments made, contact Thomas via email
directly. If these developments are within the focus of the developers'
research, then a collaboration, i.e. resulting in a co-authored peer reviewed
publication, can be an option. Otherwise, you can contract us for
paid-developments.

.. _PDF: http://msnoise.org/doc/MSNoise.pdf
