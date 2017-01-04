MSNoise - Documentation
=======================

Originally, MSNoise was a "Python Package for Monitoring Seismic Velocity
Changes using Ambient Seismic Noise". With the release of MSNoise 1.4, and
because of the Plugin Support, we could call MSNoise: "Measuring with Seismic
Noise".

The standard MSNoise workflow is designed to go from seismic data archives to
dv/v curves. The monitoring is achieved by computing the cross-correlation of
continuous seismic records for each pair of a network and by studying the
changes in the crosscorrelation function relative to a reference.

The goal of the "suite" is to provide researchers with an efficient processing
tool, while keeping the need for coding to a minimum and avoiding being a black
box. Moreover, as long as the in- and outputs of each step are respected, they
can easily be replaced with one's own codes ! (See :ref:`workflow`).

Plugins can be added and extend the standard workflow from any steps, e.g.
using MSNoise as a cross-correlation toolbox until the `stack` step, and then
branching to the workflow provided by one's plugin.


If you use MSNoise for your research and prepare publications, please consider
citing MSNoise:
**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.


Installation
------------

.. toctree::
   :maxdepth: 1

   installation

Workflow
--------

.. toctree::
   :glob:
   :maxdepth: 1

   workflow/*

Plotting
--------

.. toctree::
   :maxdepth: 2

   plotting


Development & Miscellaneous
---------------------------

.. toctree::
   :maxdepth: 2

   plugins
   api
   core
   table_def
   how_tos
   about_db_performances
   references
   contributors

.. toctree::
    :maxdepth: 1

    clickhelp/msnoise




Release Notes:

.. toctree::
   :maxdepth: 1

   releasenotes/msnoise-1.5
   releasenotes/msnoise-1.4
   releasenotes/msnoise-1.3.1
   releasenotes/msnoise-1.3
   releasenotes/msnoise-1.2.5
   releasenotes/msnoise-1.2.4
   releasenotes/msnoise-1.2.3
   releasenotes/msnoise-1.2.2
   releasenotes/msnoise-1.2.1
   releasenotes/msnoise-1.2
   releasenotes/msnoise-1.0


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
