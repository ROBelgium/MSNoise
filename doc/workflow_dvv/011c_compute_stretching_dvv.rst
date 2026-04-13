.. include:: ../configs.hrst

Compute dv/v aggregate from Stretching
---------------------------------------

This step aggregates the per-pair stretching results into network-level dv/v
statistics (mean, std, weighted mean, percentiles, etc.) using the same
aggregation engine as the MWCS dv/v step.

To run:

.. code-block:: sh

    $ msnoise cc dtt compute_stretching_dvv

The output files are written under:

.. code-block:: text

    <output_folder>/<lineage>/stretching_dvv_1/_output/<mov_stack>/dvv_<pair_type>_<comp>.nc

.. seealso::

   :doc:`009_compute_dvv` — full description of the aggregation algorithm,
   configuration parameters, and output format (shared with all three dvv methods).

.. seealso::

   **Reading these results in Python** — use :class:`MSNoiseResult <msnoise.results.MSNoiseResult>`:

   .. code-block:: python

      from msnoise.results import MSNoiseResult
      from msnoise.core.db import connect
      db = connect()
      r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                 stack=1, refstack=1,
                                 stretching=1, stretching_dvv=1)
      ds = r.get_dvv(pair_type="CC", components="ZZ", mov_stack=("1D", "1D"))

   See :ref:`msnoise_result` for the full guide and all available methods.
