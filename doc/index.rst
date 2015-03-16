MSNoise - Documentation
===============================================================================================

MSNoise is a Python Package for Monitoring Seismic Velocity Changes using Ambient Seismic Noise.
The monitoring is achieved by computing the cross-correlation of continuous seismic records for each
pair of a network and by studying the changes in the crosscorrelation function relative to a reference.

The goal of the "suite" is to provide researchers with an efficient processing tool, while keeping the
need for coding to a minimum and avoiding being a black box.
Moreover, as long as the in- and outputs of each step are respected, they can easily
be replaced with one's own codes ! (See :ref:`workflow`).

Contents:

.. toctree::
   :maxdepth: 2

   installation

   workflow
   plotting
   api
   core
   table_def
   how_tos
   about_db_performances
   references
   todos


Release Notes:

.. toctree::
   :maxdepth: 1

   releasenotes/msnoise-1.3
   releasenotes/msnoise-1.2.5
   releasenotes/msnoise-1.2.4


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
