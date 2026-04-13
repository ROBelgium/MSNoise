.. include:: ../configs.hrst

Compute WCT (Wavelet Coherence Transform)
-----------------------------------------

.. automodule::  msnoise.s08_compute_wct

.. seealso::

   **Reading these results in Python** — use :class:`MSNoiseResult <msnoise.results.MSNoiseResult>`:

   .. code-block:: python

      from msnoise.results import MSNoiseResult
      from msnoise.core.db import connect
      db = connect()
      r = MSNoiseResult.from_ids(db, ...)  # include the steps you need
      # then call r.get_wct(...)

   See :ref:`msnoise_result` for the full guide and all available methods.
