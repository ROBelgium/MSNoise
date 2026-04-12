.. include:: ../configs.hrst

MSNoise 2.0
===========

Release type: **major**

MSNoise 2.0 is a complete architectural overhaul.  The processing pipeline
is now driven by a **workflow DAG** of config sets, steps, and lineages,
replacing the 1.x single-branch flat-file approach.

.. contents::
    :local:


Config Sets and Workflow Lineages
----------------------------------

The central innovation in 2.0 is the **config set** system.  Every processing
category (``cc``, ``filter``, ``mwcs``, ``stack``, …) can have multiple
parameterised instances running in parallel.  Each unique path through the DAG
is a *lineage* — encoded as a ``/``-separated string of step names — and
becomes the folder hierarchy under ``OUTPUT/``.

See :ref:`workflow_concepts` for a full explanation.

Key consequences:

- **Multiple filter bands** in one project: create a second ``filter`` config
  set.  Both branches are processed and stored independently.
- **Multiple MWCS parameters**: create a second ``mwcs`` config set.
- **Output paths**: ``OUTPUT/preprocess_1/cc_1/filter_1/stack_1/mwcs_1/_output/…``
  — the lineage IS the path.


New CLI commands
-----------------

Config management:

.. code-block:: sh

    msnoise config list_sets
    msnoise config create_set mwcs
    msnoise config copy_set mwcs 1 mwcs 2
    msnoise config delete_set mwcs 2
    msnoise config show_set mwcs 1
    msnoise config get cc.cc_sampling_rate
    msnoise config set mwcs.2.freqmin 1.0
    msnoise config reset mwcs.2.freqmin

Station management:

.. code-block:: sh

    msnoise utils import-stationxml inventory.xml
    msnoise utils import-stationxml https://fdsn.org/.../query

Job utilities:

.. code-block:: sh

    msnoise utils create_preprocess_jobs --date 2024-01-01
    msnoise utils create_psd_jobs --date_range 2024-01-01 2024-01-31


PSD pipeline (new)
-------------------

MSNoise 2.0 adds a dedicated PSD quality-control pipeline, independent
of the CC pipeline:

- :mod:`msnoise.s20_psd_compute` — per-day PSD using ObsPy PPSD
- ``msnoise.s21_psd_compute_rms`` — per-frequency-band RMS


WCT fused mode
---------------

The Wavelet Coherence Transform step (:mod:`msnoise.s08_compute_wct`) now
supports a **fused mode** (``|wavelet.wct_compute_dtt|`` = ``True``, default):
dt/t results are computed and written directly inline, skipping intermediate
WCT file storage (~2000× storage reduction per pair/year).  The standalone
``msnoise.s09_compute_wct_dtt`` step remains available for reprocessing
existing WCT files.


DataSource abstraction
-----------------------

Stations can now be associated with individual data sources (local SDS archive
or FDSN/EIDA remote service) via the ``DataSource`` table.  This enables
mixed local+remote network processing in a single project.


Plugin API
-----------

The plugin system is now fully topology-aware.  Plugins can:

- Register new workflow categories via ``msnoise.plugins.workflow_chains``
  and ``msnoise.plugins.workflow_order`` entry points.  New categories are
  automatically valid in ``msnoise new_jobs --after``, in job propagation,
  and in ``MSNoiseResult``.
- Import from the stable :mod:`msnoise.plugins` surface instead of
  ``msnoise.core.*``.

See :doc:`../plugins` for the full guide.


Removed / renamed
------------------

- ``msnoise/move2obspy.py`` → ``msnoise/core/compute.py``
- ``msnoise/preprocessing.py`` → ``msnoise/core/preprocessing.py``
- ``msnoise/default.py`` → per-category CSV files in ``msnoise/config/``
- ``msnoise/psd_compute.py`` → ``msnoise/s20_psd_compute.py``
- ``msnoise/psd_compute_rms.py`` → ``msnoise/s21_psd_compute_rms.py``
- ``LayeredParams`` class renamed to ``MSNoiseParams``
