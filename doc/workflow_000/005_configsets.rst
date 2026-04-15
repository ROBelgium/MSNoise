.. include:: ../configs.hrst

.. _workflow_configsets:

Setting up a Multi-Branch Project
==================================

This page walks through configuring a project that uses **multiple config sets**
to run parallel processing branches — for example, two filter frequency bands
or two MWCS window lengths.

For the conceptual background, see :ref:`workflow_concepts`.

.. contents::
    :local:


Default project (one branch)
------------------------------

A freshly initialised project has one config set per category:
``cc_1``, ``filter_1``, ``stack_1``, ``refstack_1``, ``mwcs_1``, etc.
The whole pipeline runs as a single linear chain.


Adding a second filter band
-----------------------------

1. Create the second config set and edit its parameters:

   .. code-block:: sh

       msnoise config create_set filter
       msnoise config list_sets                 # should show filter_1 and filter_2
       msnoise config set filter.2.freqmin 1.0
       msnoise config set filter.2.freqmax 5.0

2. Apply the new topology to the database:

   .. code-block:: sh

       msnoise db upgrade

   This creates ``filter_2`` as a new WorkflowStep and wires it:
   ``cc_1 → filter_2 → stack_1`` and ``cc_1 → filter_2 → refstack_1`` (siblings)

3. Run the CC step as usual.  Both filter branches are processed:

   .. code-block:: sh

       msnoise cc compute
       msnoise cc stack
       msnoise cc stack_refstack
       msnoise cc dtt compute_mwcs
       …

   Output files for ``filter_1`` land under::

       OUTPUT/preprocess_1/cc_1/filter_1/stack_1/…

   Output files for ``filter_2`` land under::

       OUTPUT/preprocess_1/cc_1/filter_2/stack_1/…


Adding a second MWCS parameterisation
---------------------------------------

.. code-block:: sh

    msnoise config create_set mwcs
    msnoise config copy_set mwcs 1 mwcs 2    # start from set 1's values
    msnoise config set mwcs.2.freqmin 1.0
    msnoise config set mwcs.2.freqmax 5.0
    msnoise config set mwcs.2.mwcs_wlen 5
    msnoise db upgrade

Both ``mwcs_1`` and ``mwcs_2`` will run on every refstack branch.


Checking what was created
--------------------------

.. code-block:: sh

    msnoise config list_sets          # all categories and set numbers
    msnoise admin                     # web UI → Workflow → Steps

From Python:

.. code-block:: python

    from msnoise.plugins import connect, get_workflow_steps
    db = connect()
    for s in get_workflow_steps(db):
        print(f"{s.step_name:25s}  cat={s.category}  set={s.set_number}")


Resetting jobs after a config change
--------------------------------------

If you add a new config set after data has already been processed, reset
the relevant steps so the new branch gets processed:

.. code-block:: sh

    msnoise reset stack_1 --all      # reset all stack jobs to T
    msnoise new_jobs --after cc      # re-propagate from cc → stack → …

