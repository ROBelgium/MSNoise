.. include:: ../configs.hrst

Compute MWCS (Moving Window Cross-Spectrum)
-------------------------------------------

.. automodule::  msnoise.s05_compute_mwcs

.. seealso::

   **Reading these results in Python** — use :class:`MSNoiseResult <msnoise.results.MSNoiseResult>`:

   .. code-block:: python

      from msnoise.results import MSNoiseResult
      from msnoise.core.db import connect
      db = connect()
      r = MSNoiseResult.from_ids(db, ...)  # include the steps you need
      # then call r.get_mwcs(...)

   See :ref:`msnoise_result` for the full guide and all available methods.
