.. include:: ../configs.hrst


.. _workflowcomputedtt:

Compute dt/t from MWCS measurements
-----------------------------------

.. automodule::  msnoise.s06_compute_mwcs_dtt

.. seealso::

   **Reading these results in Python** — use :class:`MSNoiseResult <msnoise.results.MSNoiseResult>`:

   .. code-block:: python

      from msnoise.results import MSNoiseResult
      from msnoise.core.db import connect
      db = connect()
      r = MSNoiseResult.from_ids(db, ...)  # include the steps you need
      # then call r.get_mwcs_dtt(...)

   See :ref:`msnoise_result` for the full guide and all available methods.
