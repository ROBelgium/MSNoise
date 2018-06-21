.. include:: ../configs.hrst

MSNoise 1.6
===========

Release date: XX June 2018


Release type: major


Release notes:

.. contents::
    :local:

Introduction
------------
More than 1 year after the last major release (:doc:`msnoise-1.5`) I'm proud
to announce the new :doc:`msnoise-1.6`. It is a **major** release, with a
massive amount of work since the last one: in `GitHub numbers
<https://github.com/ROBelgium/MSNoise/graphs/contributors?from=2017-04-28&to=2017-06-25&type=c>`_
, it's over TODO commits and over TODO lines of code and documentation changed
or added!


MSNoise 1.6 introduces a series of **new features** :

* The workflow has been rewritten to create "job types" per step, making it 
  easier for users to reset a jobs before a specific step. 
* This and other smaller adaptations to the code allows to run MSNoise more
  effiently, e.g. on a HPC.



This version has benefited from outputs/ideas/pull requests/questions from
several users/friends (listed alphabetically):

# TODO

* Raphael De Plaen
* Clare Donaldson
* Robert Green
* Aurelien Mordret
* Lukas Preiswerk
* The participants to the NERC MSNoise Liverpool Workshop in January 2017
* all others (don't be mad :-) )


Thanks to all for using MSNoise, and please, let us know why/how you use it
(and please cite it!)!

To date, we found/are aware of 48 publications using MSNoise ! That's the best
validation of our project ever and it has doubled since last year !! 


*Thomas and friends...*


~~~~

PS: if you use MSNoise for your research and prepare publications, **please
consider citing it**:

**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.


Requirements
------------

There were no changes in the requirements. Note that MSNoise is always tested
against the latest release versions of the main packages, so older installations
that are not maintainted/updated regularly (years) could encounter issues.


Web-based Admin Interface Changes
---------------------------------

* Feature: The pagination size (100, 200, etc. row) is now allowed on the
  Station, Job and DataAvailability views are now 


Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* ADDED: ``hpc`` for flagging if MSNoise is running on HPC. Not yet used.


Top level DB command
--------------------

I've added a new command group called `db` that gathers all db-related actions:

* ``msnoise db upgrade`` is a replacement for the ``msnoise upgrade_db``
* ``msnoise db clean_duplicates`` deletes duplicate jobs (might happen). Unique
  sets of ``day``, ``pair`` and ``jobtypes``are considered.

In the future, this command will include the possibility to dump and load 
entire databases, for example for switching from a local sqlite instance to a 
large MySQL server.

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
stacks if needed, then will reset the ``STACK`` jobs to ``Todo`` so that "mov
stacks can be done. Note that calling ``msnoise stack -r -m`` will execute 
the steps in the right order.


Expert/Lazy mode
~~~~~~~~~~~~~~~~

* Added the possibility to walk in subfolders recursively by using 
  ``scan_archive --path -r``


Preprocessing and Cross-Correlation
-----------------------------------

Preprocessing
~~~~~~~~~~~~~

* Only small changes were done for this step, mainly checks of matching
  sampling rates, empty streams. DB-related optimisations make this step 
  faster too.

Cross-Correlation TODO
~~~~~~~~~~~~~~~~~

* It is now possible to do the Cross-Correlation (classic "CC"), the Auto-
  Correlation ("AC") or the Cross-Components within the same station ("SC").
  To achieve this, we removed the `ZZ`, `ZT`, etc parameters from the
  configuration and replaced it with ``components_to_compute`` which takes a list:
  e.g. `ZZ,ZE,ZN,EZ,EE,EN,NZ,NE,NN` for the full non-rotated tensor between
  stations. If `autocorr` is set to "Y", then the cross-components (SC) of each
  station will also be computed (of course, `ZE` and `EZ` are identical).

* The cross-correlation is done on sliding windows on the available data. If one
  trace contains a gap, the window is skipped. This corrects previous errors
  linked with gaps synchronised in time that lead to perfect sinc autocorr
  functions. The windows should have a duration of at least 2 times the `maxlag`
  configuration value.

* The whitening procedure can be skipped by setting the ``whitening``
  configuration to `None`. The two other ``whitening`` modes are "[A]ll except
  for auto-correlation" or "Only if [C]omponents are different". This allows
  skipping the whitening when, for example, computing ZZ components for very
  close by stations (much closer than the wavelength sampled), leading to
  spatial autocorrelation issues.


Command Line changes
--------------------


* ``msnoise db`` a new top-level group command for gathering the db-related 
  commands.
* ``msnoise db clean_duplicates`` allows to remove duplicates from the jobs 
  table.

* ``msnoise scan_archive --path -r`` will only scan the ``path`` 
independently of its structure. Passing the ``-r`` it will walk in subfolders.

Note, all commands are documented: :doc:`../clickhelp/msnoise`.


API Changes
-----------

* New ``get_params`` API method to return a ``Params`` class containing all 
  configuration bits, avoiding unnecessary calls to the DB later.
* Many (many) small optimizations to the core functions to make them faster,
  mostly when interacting with the DB, thanks to ``get_params``.
* ``update_config`` can now modify parameters for installed plugins.
* New ``massive_update_job`` to update a list of ``Jobs`` to a given flag.
* ``build_ref_datelist`` and ``build_movstack_datelist`` return the smallest 
  date from the data_availability table if the ``ref_begin`` or ``startdate``
  have not been modified from ``1970-01-01``.
* Removed ``linear_regression`` as this is now included in ObsPy.
* Modified ``get_dtt_next_job`` returns jobs in random order (will maybe change
  again later when the HPC mode gets more developed.
* Added missing documentation for several methods.

See :doc:`../api`.


Other Bugfixes TODO
--------------

* Removed the call to ``scipy.stats.nanmean`` and replaced by ``numpy.nanmean``
* Better error message in compute_cc when the content of the slice is only zeros
  or smaller than ``rms_threshold``
* Checked all "integer" - "float" warnings from numpy/scipy
* crondays were hardcoded to -2, now taking the ``crondays`` value from the DB
* Py3 error when ``msnoise scan_archive`` in cron mode


Plot Updates TODO
------------

* ``msnoise plot ccftime`` now accepts -e (--envelope) and will plot the
  envelope of the ccfs.
* ``msnoise plot ccftime``, ``msnoise plot interferogram`` and
  ``msnoise plot distance`` now accept -r (--refilter) to refilter the CCFs
  before plotting. The argument must be a column-separated string (e.g.
  ``-r 4:8`` for filtering between 4.0 and 8.0 Hz).
* ``msnoise plot distance`` accepts a new ``--virtual-source`` NET.STA parameter
  to only plot the pairs including this station.
* Most plots have better titles (filter details, etc).
* The dv/v plot now allows averaging over components by passing them as comma-
  separated values.

See :doc:`../plotting`.

Performance and Code improvements TODO
---------------------------------

Improvements in terms of performances have also been done for MSNoise 1.5:

* Added fftpack optimized nfft (scipy's next_fast_len). This could lead to some
  small differences in the final result of the MWCS procedure, because of the
  number of points used for smoothing the (cross-)spectra.
* Replaced binarization (sign) and windsorizing (clip) by standard numpy
  functions operating directly inplace on the arrays, avoiding unecessary
  copies.
* The preprocessing only reads files that should contain the right component.
* The stretching code has been improved (thanks to Clare Donaldson)

Documentation TODO
-------------

* New tutorial for installing and configuring MySQL and MySQL workbench.
* The Workflow is separated in steps (See :ref:`workflow`))
* The Documentation is now available in PDF format (PDF_).

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