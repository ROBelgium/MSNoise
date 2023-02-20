MSNoise - Documentation
=======================

Originally, MSNoise was a "Python Package for Monitoring Seismic Velocity
Changes using Ambient Seismic Noise". With the release of MSNoise 1.4, and
because of the Plugin Support, we could call MSNoise: "Measuring with Seismic
Noise".

The current release version of MSNoise is **MSNoise 2.0** and is a
Python 3 version only.

The standard MSNoise workflow is designed to go from seismic data archives to
dv/v curves. The monitoring is achieved by computing the cross-correlation of
continuous seismic records for each pair of a network and by studying the
changes in the cross-correlation function relative to a reference.

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

This documentation is also available in PDF format on the MSNoise Website
(PDF_).

Installation
============

.. toctree::
    :maxdepth: 2

    installation

Workflows:
==========

Initializing and Configuring a Project
--------------------------------------
Commands in this workflow have the ``msnoise <command>`` form.

.. toctree::
    :maxdepth: 2

    workflow_000/0000_intro.rst


Cross-Correlation Functions
---------------------------
Commands in this workflow have the ``msnoise cc <command>`` form.

.. toctree::
    :maxdepth: 2

    workflow_ccf/0000_intro.rst

Computing Relative Velocity Variations
--------------------------------------
Commands in this workflow have the ``msnoise cc dvv <command>`` form.

.. toctree::
    :maxdepth: 2

    workflow_dvv/0000_intro.rst

Computing Power Spectral Densities (Quality Control)
----------------------------------------------------
Commands in this workflow have the ``msnoise qc <command>`` form.

.. toctree::
    :maxdepth: 2

    workflow_psd/0000_intro.rst


Interacting with MSNoise
========================


.. toctree::
    :maxdepth: 2
    :glob:

    how_tos
    auto_examples/index
    api
    core
    plugins
    

.. toctree::
    :maxdepth: 1

    clickhelp/msnoise

Development & Miscellaneous
===========================

.. toctree::
    :maxdepth: 2

    table_def
    about_db_performances
    references
    contributors
    releasenotes


.. _PDF: http://msnoise.org/doc/MSNoise.pdf