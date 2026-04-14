MSNoise - Documentation
=======================

MSNoise is a Python package for monitoring seismic velocity changes using
ambient seismic noise.  It automates the full pipeline from raw seismic
archives to dv/v curves: cross-correlations, stacking, reference computation,
and relative velocity estimation via MWCS, stretching, or wavelet methods.

MSNoise 2.x is a full rewrite built on xarray, with a config-set / lineage
architecture that makes it easy to run the same data through multiple
parameter branches without reprocessing from scratch.  Plugins can extend
or replace any step.

If you use MSNoise for your research, please cite:

| **Lecocq, T., C. Caudron, et F. Brenguier (2014)**, MSNoise, a Python Package
| for Monitoring Seismic Velocity Changes Using Ambient Seismic Noise,
| *Seismological Research Letters*, 85(3), 715–726, doi:10.1785/0220130073.


.. toctree::
    :maxdepth: 2
    :caption: Getting Started

    installation
    workflow_concepts

Install MSNoise and learn the 2.x concepts — config sets, lineages, and
jobs — before running anything.


.. toctree::
    :maxdepth: 2
    :caption: Workflows

    workflow_000/0000_intro.rst
    workflow_ccf/0000_intro.rst
    workflow_dvv/0000_intro.rst
    workflow_psd/0000_intro.rst

Step-by-step guides from project initialisation to dv/v and PSD outputs,
with the exact commands at each stage.


.. toctree::
    :maxdepth: 2
    :caption: Using MSNoise

    results
    how_tos
    auto_examples/index

How to read and work with computed results in notebooks and scripts, plus
recipes for common tasks and a gallery of worked examples.


.. toctree::
    :maxdepth: 2
    :caption: API Reference

    core
    api
    table_def
    plugins

Full reference for the Python API, internal modules, database schema, and
the plugin system.


.. toctree::
    :maxdepth: 1
    :caption: CLI Reference

    clickhelp/msnoise

Every ``msnoise`` command, its arguments, and options.


.. toctree::
    :maxdepth: 2
    :caption: Development

    about_db_performances
    references
    contributors
    releasenotes


.. _PDF: http://msnoise.org/doc/MSNoise.pdf