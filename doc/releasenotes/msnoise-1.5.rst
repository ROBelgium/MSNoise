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
Just over XXX after the last major release (:doc:`msnoise-1.4`) we are proud
to announce the new :doc:`msnoise-1.5`. It is a **major** release, with a
massive amount of work since the last one: in `GitHub numbers
<https://github.com/ROBelgium/MSNoise/graphs/contributors?from=2016-06-02&to=2017-04-30&type=c>`_
, it's over 100 commits and about 2500 lines of code and documentation changed
or added!

MSNoise 1.5 introduces **XXXX major new features** :

* new preprocessing function, returns streams, NO MERGE TO ZEROS anymore
* started to move core math functions to obspy (and move2obspy)
* default Lanczos resampling method now
* all-component possible
* autocorr-not-whitened + more explanation from lukas about whitening types
* dvv plot allows averaging over components
* stretching improvements
* reduced the memory imprint caused by scipy's caching of FFT arrays
* cross-correlation between slices are only possible when the two traces are longer than 2 times the maxlag defined
* add a default delay of 1 (customizable) to start parallel threads (using -t)
* new `make_same_length` API method to return common data for two traces
* this documentation is now available in PDF too (easier for offline usage)
* custom archive structure can now directly be input in the config (need at least one slash)
* everything in compute_cc and preprocessing is obspy stream
*


Parameters

* REMOVED: decimation_factor : it's now computed automatically
* REMOVED: ZZ, RR, TT, etc in favour of `components_to_compute`
* ADDED: components_to_compute: it's accepting a comma-separated list of components
* ADDED: whitening method +LUKAS


In MSNoise 1.5 we drop the dependency to `statsmodels` in favour of our own
linear regression routine that we included in `obspy` as of version 1.1.


This version has benefited from outputs/ideas/pull requests/questions from
several users/friends:

* Lukas Preiswerk
* Clare Donaldson
* Robert Green
* ...
* all others (don't be mad :-) )


Thanks to all for using MSNoise, and please, let us know why/how you use it
(and please cite it!)!

To date, we found/are aware of XXX publications using MSNoise ! That's the best
validation of our project ever ! See the full list on the
`MSNoise website <http://www.msnoise.org/they-cite-msnoise/>`_.


*Thomas & Corentin*


~~~~

PS: if you use MSNoise for your research and prepare publications, please
consider citing it:

**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.


Requirements
------------

If you have any package older than those mentioned here, some things might not work as expected.

* Pandas > 0.18.0
* obspy > 1.1
* Scikits-samplerate no longer needed
* statsmodels is no longer needed


Web-based Admin Interface Changes
---------------------------------

* Bugfix: Stations can be added manually
* Bugfix: Bulk operations on Job based on their jobtype is now possible
* Feature: Calling admin/data_availability.png should serve the image directy
* Bugfix: Set debug-mode to False to increase performance and remove the risk
  of Werkzeug failing.
* Feature: It's now possible to "Brand" MSNoise by setting an environment
  variable `msnoise_brand` to `CustomName|http://custom.url/logo.png`



Plugin support
--------------

* Plugins can now declare their own templates using the
  `msnoise.plugins.templates` entry point.
* Plugins can override the "components_to_compute" config bit and the MSNoise
  API `get_components_to_compute(session, plugin=None)` works like `get_config`.

See :doc:`../plugins`.


Command Line changes
--------------------


* msnoise config --sync
* msnoise scan_archive --path
* msnoise populate --fromDA
* msnoise p == msnoise plugin
* add a default delay of 1 (customizable) to start parallel threads (using -t)


Note, all commands are documented: :doc:`../clickhelp/msnoise`.


Other Bugfixes
--------------

* Removed the call to `scipy.stats.nanmean` and replaced by `numpy.nanmean`
* Better error message in compute_cc when the content of the slice is only zeros
  or smaller than `rms_threshold`
* Checked all "integer" - "float" warnings from numpy/scipy
* crondays were hardcoded to -2, now taking the `crondays` value from the DB
* Py3 error when `msnoise scan_archive` in cron mode




Plot Updates
------------

* The ccftime now accepts -e (--envelope) and will plot the envelope of the ccfs
* Most plots have better titles (filter details, etc)


Performance improvements
------------------------

Improvements in terms of performances have also been done for MSNoise 1.5:

* Added fftpack optimized nfft (scipy's next_fast_len) !! smoothing warning !!
* Replaced binarization (sign) and windsorizing (clip) by standard numpy functions
  operating directly inplace on the arrays, avoiding unecessary copies
* preprocessing only reads files that should contain the right component



Upgrading an existing project to MSNoise 1.5
--------------------------------------------

Some users will want to keep their current project without recomputing
everything. This requires adding a few configuration parameters to the database

Running the following command will take care of the upgrade from 1.4 to 1.5:

.. code-block:: sh

    msnoise upgrade_db





A final note about development pace and choices
-----------------------------------------------

* MSNoise team is

  * **1 developper** (Thomas)
  * 1 dedicated debugger (Corentin)
  * less than 5 really *active* users, providing feedback and/or lines of codes
    (Esteban, Raphaël, Aurélien, Carmelo, ...)

* All software engineering ideas are coming from too infrequent beerstormings
  between Thomas & others
* The web-interface and the plugin support were developed during Thomas'
  holidays

If you need help, please ask your question on the mailing list. Don't be afraid
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
