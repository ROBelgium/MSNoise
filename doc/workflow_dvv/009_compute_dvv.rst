.. include:: ../configs.hrst

Compute dv/v from dt/t(MWCS)
----------------------------

.. automodule::  msnoise.s07_compute_dvv

.. seealso::

   **Reading these results in Python** — use :class:`MSNoiseResult <msnoise.results.MSNoiseResult>`:

   .. code-block:: python

      from msnoise.results import MSNoiseResult
      from msnoise.core.db import connect
      db = connect()
      r = MSNoiseResult.from_ids(db, ...)  # include the steps you need
      # then call r.get_dvv(...)

   See :ref:`msnoise_result` for the full guide and all available methods.
