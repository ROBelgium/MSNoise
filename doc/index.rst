MSNoise - Documentation
=======================

MSNoise is a Python Package for Monitoring Seismic Velocity Changes using Ambient Seismic Noise.
The monitoring is achieved by computing the cross-correlation of continuous seismic records for each
pair of a network and by studying the changes in the crosscorrelation function relative to a reference.

The goal of the "suite" is to provide researchers with an efficient processing tool, while keeping the
need for coding to a minimum and avoiding being a black box.
Moreover, as long as the in- and outputs of each step are respected, they can easily
be replaced with one's own codes ! (See :ref:`workflow`).

If you use MSNoise for your research and prepare publications, please consider
citing MSNoise:
**Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
*Seismological Research Letters*, 85(3), 715‑726, doi:10.1785/0220130073.


Documentation Content:

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

....

Release Notes:

.. toctree::
   :maxdepth: 1

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
