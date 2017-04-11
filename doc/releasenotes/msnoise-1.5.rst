.. include:: ../configs.hrst

MSNoise 1.5
===========

Release date: XX XXXXX XXXX


Release type: major


Release notes:

.. contents::
    :local:

Introduction
------------
About 1 year after the last major release (:doc:`msnoise-1.4`) we are proud
to announce the new :doc:`msnoise-1.5`. It is a **major** release, with a
massive amount of work since the last one: in `GitHub numbers
<https://github.com/ROBelgium/MSNoise/graphs/contributors?from=2016-06-02&to=2017-04-30&type=c>`_
, it's over 140 commits and about 2500 lines of code and documentation changed
or added!


MSNoise 1.5 introduces a series of **new features** :

* We have started to move core math functions to ObsPy, currently the only one
  ready is `linear_regression`, a function I wrote to remove the dependency to
  ``statsmodels``, required to move `mwcs` to ObsPy later.

* The preprocessing routine has been isolated, rewritten and optimized. It is
  now a standalone script, callable by plugins. It returns a Stream object with
  all the data needed for the analysis.

* This change in preprocessing was done mostly to allow cross-component, auto-
  correlation and cross-correlation, with or without rotation, to be done with
  the same code. CC, SC and AC are now supported in MSNoise with proper
  whitening (possible to disable spectral whitening for specific cases).

* This documentation is now available in PDF_ too (easier for offline usage) and
  it also includes a new tutorial for setting up the MySQL server and Workbench.

* Last but not least: MSNoise is "tested" automatically on Linux (thanks to
  TravisCI) & Windows (thanks to Appveyor), for Python versions 2.7 and 3.5.
  With MSNoise 1.5 we also added the MacOSX tests on TravisCI. With these tests,
  we can guarantee MSNoise works on different platforms and Anaconda
  (or miniconda) python versions.


This version has benefited from outputs/ideas/pull requests/questions from
several users/friends (listed alphabetically):

* Raphael De Plaen
* Clare Donaldson
* Robert Green
* Aurelien Mordret
* Lukas Preiswerk
* The participants to the NERC MSNoise Liverpool Workshop in January 2017
* all others (don't be mad :-) )


Thanks to all for using MSNoise, and please, let us know why/how you use it
(and please cite it!)!

To date, we found/are aware of 25 publications using MSNoise ! That's the best
validation of our project ever ! See the full list on the
`MSNoise website <http://www.msnoise.org/they-cite-msnoise/>`_.


*Thomas, Corentin and others*


~~~~

PS: if you use MSNoise for your research and prepare publications, **please
consider citing it**:

**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.


Requirements
------------

If you have any package older than those mentioned here, some things might not
work as expected.

* Pandas >= 0.18.0
* obspy >= 1.1

No longer needed:

* scikits-samplerate
* statsmodels


Web-based Admin Interface Changes
---------------------------------

* Bugfix: Stations can now be added manually
* Bugfix: Bulk operations on Job based on their jobtype is now possible
* Bugfix: Set debug-mode to False to increase performance and remove the risk
  of Werkzeug failing.

* Feature: Calling admin/data_availability.png should serve the image directly
* Feature: It's now possible to "Brand" MSNoise:
  The name and the logo of the page can be overriden by setting an environment
  variable with a name and the HTML tag of the logo image:

  .. code:: sh

      set msnoise_brand="ROB|<img src='http://www.seismologie.be/img/oma/ROB-logo.svg' width=200 height=200>"

  and then starting msnoise admin:

  .. image:: ../.static/branding.png
      :align: center

Configuration Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

* REMOVED: ``decimation_factor``: it's now computed automatically from the
  stream's sampling rate. If the decimation factor is non-integer, the routine
  exits and the user is advised to use the Lanczos resampler.
* REMOVED: `ZZ`, `RR`, `TT`, etc in favor of ``components_to_compute``.
* ADDED: ``components_to_compute``: it's accepting a comma-separated list of
  components.
* ADDED: ``whitening`` method for disabling spectral whitening in specific
  cases.

Populate Station Table and Scan Archive
---------------------------------------


* Because of platform-dependent issues, we dropped the need for the "find"
  command to be present on each machine. The "cron" scan_archive procedure is
  now a simple python routine. It might take a little longer, but much safer.
  This also means the ``crondays`` parameter can be set to a floating point
  value, "-1.0" meaning "files modified in the last 1 day".

* The custom archive structure can now directly be input in the config.
  This is still quite preliminary and needs at least one slash in the parameter
  value to be identified as such by the code.

Expert/Lazy mode
~~~~~~~~~~~~~~~~

* Added an Expert/Lazy mode to scan files directly by passing the path of their
  containing folder to ``scan_archive --path`` :ref:`(see here)<scan-archive-expert>`.

* Added an Expert/Lazy mode to input network/stations directly from the data
  scanned using the new ``populate --fromDA`` procedure
  (:ref:`(see here)<populate-expert>`).

* The ``scan_archive`` is now capable of handling multiplexed files and inputs
  one line in DataAvailability per unique "trace id" in files.

Practically, the three changes above allow users to:

#. Download a single or few mseed file with N stations, for example from webdc
   or using FDSN webservices in ObsPy or else.
#. Run ``msnoise install``
#. Run ``msnoise scan_archive --path /path/to/the/big/file/``
#. Run ``msnoise populate --fromDA``
#. Run ``msnoise new_jobs --init``
#. Run ``msnoise admin``, define (pre-)processing and filter parameters
#. Run ``msnoise compute_cc``


New Jobs
--------

* The ``msnoise new_jobs`` routine has been rewritten to correct a few bugs
  linked to start/end dates (missed days if a file of day 2 begins on day 1 and
  ends on day 3, jobs were only created for day 1 and 3. This is corrected.

* A change in the database structure (added an index on the
  "day"+"pair"+"jobtype") allows running ``msnoise new_jobs`` much faster than
  before. It took 5m40s for a test with 105381 jobs, while it would have taken
  hours before. For the first run, it is still faster (20s for the example
  above) to use ``msnoise new_jobs --init`` as this doesn't check for existing
  jobs.

.. warning::

    The change in the database (adding an index) requires that you
    ``msnoise upgrade_db`` every project!


Preprocessing and Cross-Correlation
-----------------------------------

Preprocessing
~~~~~~~~~~~~~

* A complete rewrite of the preprocessing function to avoid padding and merging
  with zeros. The preprocessing function is now separated from the
  ``compute_cc`` code and can be called by external plugins. It returns a Stream
  object that can easily be filtered or slided.
* The instrument response removal has been accelerated by doing it after
  decimation, on fewer data points (Thanks to Robert Green).
  The response removal is done without the `evalresp` stuff from ObsPy, it's
  faster but potentially a little less safe.
* The default decimation tool is now Lanczos (builtin in ObsPy) and
  scikits.samplerate is no longer needed.

Cross-Correlation
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


* ``msnoise config --sync`` will try to parse the dataless files for existing
  stations in the Station table and, if found, will input the coordinates.
* ``msnoise scan_archive --path`` will only scan the ``path`` independently of
  its structure. It will only read files, not "walk" in subfolders.
* ``msnoise populate --fromDA`` will populate the station table from the
  existing data_availability table.
* ``msnoise p`` is a lazy alias to ``msnoise plugin``
* add a default delay of 1 second (customizable) to start parallel
  threads (using -t)


Note, all commands are documented: :doc:`../clickhelp/msnoise`.


API Changes
-----------

* New ``make_same_length`` API method to return common data for two traces, this
  is necessary to compute the rotation of the horizontal traces to Radial/
  Transverse.

* New ``clean_scipy_cache`` API method to reduce the memory imprint caused by
  scipy's automatic caching of FFT arrays.

See :doc:`../api`.

Plugin support
--------------

* Plugins can now declare their own templates using the
  `msnoise.plugins.templates` entry point.
* Plugins can override the "components_to_compute" config bit and the MSNoise
  API `get_components_to_compute(session, plugin=None)` works like `get_config`.

See :doc:`../plugins`.


Other Bugfixes
--------------

* Removed the call to ``scipy.stats.nanmean`` and replaced by ``numpy.nanmean``
* Better error message in compute_cc when the content of the slice is only zeros
  or smaller than ``rms_threshold``
* Checked all "integer" - "float" warnings from numpy/scipy
* crondays were hardcoded to -2, now taking the ``crondays`` value from the DB
* Py3 error when ``msnoise scan_archive`` in cron mode


Plot Updates
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

Performance and Code improvements
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

Documentation
-------------

* New tutorial for installing and configuring MySQL and MySQL workbench.
* The Workflow is separated in steps (See :ref:`workflow`))
* The Documentation is now available in PDF format (PDF_).

Upgrading an existing project to MSNoise 1.5
--------------------------------------------

Some users will want to keep their current project without recomputing
everything. This requires adding a few configuration parameters to the database

Running the following command will take care of the upgrade from 1.4 to 1.5:

.. code-block:: sh

    msnoise upgrade_db


.. warning::

    Upgradind the database **will not** remove deprecated configuration bits, so
    users should remember to define, for example, the ``components_to_compute``
    parameter if anything else than ``ZZ`` was set before.


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