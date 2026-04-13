.. include:: configs.hrst

.. _workflow_concepts:

*******************************
Workflow Concepts (MSNoise 2.x)
*******************************

MSNoise 2.x introduced a fundamentally different architecture from 1.x.
This page explains the key concepts — **config sets**, **workflow steps**,
**lineages**, and **jobs** — that underpin everything else.

.. contents::
    :local:
    :depth: 2


.. _concepts_configsets:

Config Sets
===========

A *config set* is a named collection of parameters for one processing category.
Every category (``cc``, ``mwcs``, ``filter``, ``stack``, …) can have **one or more**
config sets, numbered starting at 1.

You can think of a config set as "one specific parameterisation of one processing step".

Creating and inspecting config sets
------------------------------------

.. code-block:: sh

    # List all config sets
    msnoise config list_sets

    # Create a second mwcs config set (inherits defaults)
    msnoise config create_set mwcs

    # Inspect and edit
    msnoise config list mwcs          # all mwcs sets
    msnoise config list mwcs.2        # mwcs set 2 only
    msnoise config set mwcs.2.freqmin 1.0
    msnoise config set mwcs.2.freqmax 5.0

    # Copy an existing set as a starting point
    msnoise config copy_set mwcs 1 mwcs 2

Parameter notation uses dots: ``category.set_number.param_name``.
For set 1 (the default), the set number can be omitted: ``mwcs.freqmin``.

Practical example: two filter bands
-------------------------------------

Create a second ``filter`` config set for a higher frequency band:

.. code-block:: sh

    msnoise config create_set filter
    msnoise config set filter.2.freqmin 1.0
    msnoise config set filter.2.freqmax 5.0

After running ``msnoise db upgrade`` (or project initialisation), MSNoise creates
a second CC branch through ``filter_2``.  Both filter branches share the same
preprocessed data and CC results, but produce independent stacked CCFs,
MWCS results, and dv/v curves — all written to separate folders under ``OUTPUT/``.


.. _concepts_steps_links:

Workflow Steps and Links
=========================

Every config set becomes a **WorkflowStep** in the database (table
``WorkflowStep``).  Steps are then wired together by **WorkflowLinks**
(table ``WorkflowLink``), which form the processing DAG.

The topology follows the canonical order:

.. code-block:: text

    global_1
    ├── preprocess_1  →  cc_1  →  filter_1  →  stack_1  →  refstack_1
    │                              filter_2  →  stack_1 ↗
    │                                         ↓
    │                                 mwcs_1  →  mwcs_dtt_1  →  mwcs_dtt_dvv_1
    │                                 stretching_1            →  stretching_dvv_1
    │                                 wavelet_1  →  wavelet_dtt_1  →  wavelet_dtt_dvv_1
    └── psd_1  →  psd_rms_1

Steps are created automatically from config sets.  The admin web UI
(``msnoise admin``) lets you view and edit the workflow graph.

Pass-through categories
------------------------

Some categories (``filter``, ``global``) are *pass-through* nodes: they have
no worker script and no jobs.  They exist purely as parameter namespaces.
``propagate_downstream`` recurses through them transparently.


.. _concepts_lineages:

Lineages
=========

A *lineage* encodes the full path from the root of the DAG to the current step.
It is stored as a ``/``-separated string of step names:

.. code-block:: text

    preprocess_1/cc_1/filter_1/stack_1/refstack_1/mwcs_1

Every job carries a lineage ID that resolves to this string.  The lineage
serves two purposes:

1. **Output path**: the lineage string becomes the folder hierarchy under
   ``OUTPUT/``.  For example, MWCS results for the ``filter_1`` branch land at::

       OUTPUT/preprocess_1/cc_1/filter_1/stack_1/refstack_1/mwcs_1/_output/

   while the same step on the ``filter_2`` branch lands at::

       OUTPUT/preprocess_1/cc_1/filter_2/stack_1/refstack_1/mwcs_1/_output/

   This means two config set branches never overwrite each other.

2. **Parameter resolution**: when a worker loads its parameters, it reads the
   config for every step in its lineage, building a layered ``MSNoiseParams``
   object.  ``params.mwcs.freqmin`` therefore always refers to the *correct*
   config set for this particular branch.

Reading lineages
-----------------

.. code-block:: python

    from msnoise.plugins import connect, get_next_lineage_batch

    db = connect()
    batch = get_next_lineage_batch(db, "mwcs", group_by="pair_lineage")
    print(batch["lineage_str"])        # e.g. "preprocess_1/.../mwcs_1"
    print(batch["lineage_names"])      # ["preprocess_1", "cc_1", ..., "mwcs_1"]
    print(batch["params"].mwcs.freqmin)   # correct for this branch


.. _concepts_jobs:

Jobs
=====

A *job* is one unit of work: ``(day, pair, step, lineage)``.  Jobs are stored
in the ``Job`` table with a flag:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Flag
     - Meaning
   * - ``T``
     - **T**\odo — ready to be claimed by a worker
   * - ``I``
     - **I**\n progress — claimed atomically by one worker
   * - ``D``
     - **D**\one — completed successfully
   * - ``F``
     - **F**\ailed — worker raised an exception

Job lifecycle
--------------

.. code-block:: text

    msnoise new_jobs          →  creates T jobs for preprocess and psd
    worker claims job         →  T → I   (atomic, safe for parallel workers)
    worker finishes           →  I → D
    propagate_downstream      →  creates T jobs for the next step(s)
    worker raises exception   →  I → F

    msnoise reset cc_1        →  resets I/F jobs back to T
    msnoise reset cc_1 --all  →  resets ALL cc_1 jobs to T (including D)

HPC mode
---------

When ``|global.hpc|`` is ``Y``, ``propagate_downstream`` is **not** called
inline by workers.  Instead, the operator runs ``msnoise new_jobs --after X``
manually between steps to trigger downstream job creation.


.. _concepts_output_paths:

Output Paths
=============

All output files follow a predictable hierarchy::

    OUTPUT / lineage_names_upstream / step_name / _output / mov_stack / component / pair.nc

where ``lineage_names_upstream`` is the lineage *excluding* the current step.

Examples:

.. code-block:: text

    OUTPUT/preprocess_1/cc_1/filter_1/_output/all/ZZ/YA.UV05.00_YA.UV06.00/2024-01-01.nc
    OUTPUT/preprocess_1/cc_1/filter_1/stack_1/_output/1D_1D/ZZ/YA.UV05.00_YA.UV06.00.nc
    OUTPUT/preprocess_1/cc_1/filter_1/stack_1/refstack_1/_output/REF/ZZ/YA.UV05.00_YA.UV06.00.nc
    OUTPUT/preprocess_1/cc_1/filter_2/stack_1/refstack_1/mwcs_1/_output/1D_1D/ZZ/YA.UV05.00_YA.UV06.00.nc
    OUTPUT/psd_1/_output/daily/YA.UV05.00.HHZ/2024-01-01.nc
    OUTPUT/psd_1/psd_rms_1/_output/YA.UV05.00.HHZ.nc


.. _concepts_reading_results:

Reading Results: MSNoiseResult
==============================

Once the pipeline has run you do not need to know the on-disk path layout or
call low-level ``core.io`` functions.  The recommended interface is
:class:`MSNoiseResult <msnoise.results.MSNoiseResult>`, a single object that:

* knows which lineage branch it covers;
* exposes **only** the ``get_*`` methods valid for that branch (invalid ones
  raise ``AttributeError`` and are hidden from tab-completion);
* returns :mod:`xarray` Dataset or DataArray objects that you can slice,
  plot, or convert to pandas with one call.

.. code-block:: python

    from msnoise.results import MSNoiseResult
    from msnoise.core.db import connect

    db = connect()

    # Stacked CCFs
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
    da = r.get_ccf("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ", ("1D", "1D"))

    # Raw (pre-stack) CC outputs — only cc in lineage needed
    r_cc = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1)
    da = r_cc.get_ccf_raw("BE.UCC..HHZ:BE.MEM..HHZ", "ZZ",
                           date="2023-01-01", kind="all")

    # dv/v via MWCS
    r_dvv = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                   stack=1, refstack=1,
                                   mwcs=1, mwcs_dtt=1, mwcs_dtt_dvv=1)
    ds = r_dvv.get_dvv(pair_type="CC", components="ZZ", mov_stack=("1D", "1D"))

    # PSDs
    r_psd = MSNoiseResult.from_ids(db, psd=1, psd_rms=1)
    ds = r_psd.get_psd("BE.UCC..HHZ", day="2023-01-01")

Every workflow step page links back to the full guide.

.. seealso::

   :ref:`msnoise_result` — full ``MSNoiseResult`` guide with all methods,
   the ``kind="all"``/``"daily"`` CC output modes, branch navigation, and
   dv/v export with provenance.


Common recipes
==============

Two filter frequency bands
--------------------------

.. code-block:: sh

    msnoise config create_set filter
    msnoise config set filter.1.freqmin 0.1
    msnoise config set filter.1.freqmax 1.0
    msnoise config set filter.2.freqmin 1.0
    msnoise config set filter.2.freqmax 5.0
    msnoise db upgrade

Two MWCS window lengths
-----------------------

.. code-block:: sh

    msnoise config create_set mwcs
    msnoise config set mwcs.1.mwcs_wlen 10
    msnoise config set mwcs.2.mwcs_wlen 20
    msnoise db upgrade

Multiple moving stack windows
------------------------------

Moving stacks are configured within a single ``stack`` config set as a
tuple-of-tuples.  No additional config set is needed:

.. code-block:: sh

    msnoise config set stack.mov_stack "(('1D','1D'),('7D','1D'),('30D','1D'))"

Checking the workflow graph
----------------------------

.. code-block:: sh

    msnoise admin   # opens the web UI at http://localhost:5000
    # → Workflow → Steps / Links

or from Python:

.. code-block:: python

    from msnoise.plugins import connect, get_workflow_steps, get_workflow_chains
    db = connect()
    for step in get_workflow_steps(db):
        print(step.step_name, step.category, step.set_number)
    print(get_workflow_chains())   # full topology (plugin-aware)
