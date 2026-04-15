"""
DVV Aggregate Step
==================

Aggregates per-pair dv/v time series (from ``mwcs_dtt``, ``stretching``, or
``wavelet_dtt``) into **network-level statistics** across all station pairs.

For each ``(mov_stack, component, pair_type)`` combination, the step reads
all per-pair output files and, at each time step, computes across-pair:

- **mean** and **std** (always)
- **weighted_mean** and **weighted_std** (if ``|mwcs_dtt_dvv.dvv_weighted_mean|``
  is ``Y``; weights = 1/σ²)
- **trimmed_mean** and **trimmed_std** (if ``|mwcs_dtt_dvv.dvv_trimmed_mean|``
  is ``Y``; sigma-clip at ``|mwcs_dtt_dvv.dvv_trim_limit|``)
- **percentiles** (if ``|mwcs_dtt_dvv.dvv_percentiles|`` is set)

Pairs are split by type (CC/SC/AC) when ``|mwcs_dtt_dvv.dvv_split_pair_type|``
is ``Y``, and by component when ``|mwcs_dtt_dvv.dvv_split_components|`` is ``Y``.
Setting either to ``N`` aggregates all pairs/components into a single ``ALL`` file.

Only pairs with data quality above ``|mwcs_dtt_dvv.dvv_quality_min|`` contribute.
If ``|mwcs_dtt_dvv.dvv_output_percent|`` is ``Y``, output is in percent (×100).

To run this step:

.. code-block:: sh

    $ msnoise cc dtt compute_mwcs_dtt_dvv      # for MWCS
    $ msnoise cc dtt compute_stretching_dvv    # for Stretching
    $ msnoise cc dtt compute_wavelet_dtt_dvv   # for WCT

Aggregates per-pair dv/v results (from ``mwcs_dtt``, ``stretching``, or
``wavelet_dtt``) into network-level statistics across station pairs.

One worker handles all three DVV step categories; the category is passed as
an argument to ``main()``.  Each category maps to a parent DTT step:

====================  ===============
DVV step category     Parent category
====================  ===============
``mwcs_dtt_dvv``      ``mwcs_dtt``
``stretching_dvv``    ``stretching``
``wavelet_dtt_dvv``       ``wavelet_dtt``
====================  ===============

Output files live at::

    <root>/<lineage>/<dvv_step>/_output/<mov_stack>/dvv_<pair_type>_<comp>.nc

where ``pair_type`` is one of ``CC``, ``SC``, ``AC``, ``ALL`` and ``comp``
is the component pair string (e.g. ``ZZ``) or ``ALL``.

Configuration Parameters
------------------------

* |mwcs_dtt_dvv.dvv_split_pair_type|
* |mwcs_dtt_dvv.dvv_split_components|
* |stretching_dvv.dvv_split_pair_type|
* |stretching_dvv.dvv_split_components|
* |wavelet_dtt_dvv.dvv_split_pair_type|
* |wavelet_dtt_dvv.dvv_split_components|
* |stack.mov_stack|
* |cc.components_to_compute|
* |cc.components_to_compute_single_station|
"""
import time
import numpy as np

from .core.db import connect, get_logger
from .core.stations import get_station_pairs
from .core.workflow import get_next_lineage_batch, is_next_job_for_step, massive_update_job
from .core.io import aggregate_dvv_pairs, xr_save_dvv_agg, xr_save_dvv_pairs

# Maps DVV step category → parent DTT step category
PARENT_CATEGORY = {
    "mwcs_dtt_dvv":   "mwcs_dtt",
    "stretching_dvv": "stretching",
    "wavelet_dtt_dvv":    "wavelet_dtt",
}


def main(step_category: str = "mwcs_dtt_dvv", loglevel: str = "INFO"):
    """Aggregate per-pair dv/v into network-level statistics.

    :param step_category: One of ``"mwcs_dtt_dvv"``, ``"stretching_dvv"``,
        ``"wavelet_dtt_dvv"``.
    :param loglevel: Logging level string.
    """
    if step_category not in PARENT_CATEGORY:
        raise ValueError(
            f"Unknown DVV step category {step_category!r}. "
            f"Expected one of: {list(PARENT_CATEGORY)}"
        )
    parent_category = PARENT_CATEGORY[step_category]

    logger = get_logger(f"msnoise.{step_category}", loglevel, with_pid=True)
    logger.info(f"*** Starting: {step_category} ***")

    db = connect()

    while is_next_job_for_step(db, step_category=step_category):
        logger.debug("Getting the next batch")
        batch = get_next_lineage_batch(
            db, step_category=step_category,
            group_by="pair_lineage", loglevel=loglevel,
        )
        if batch is None:
            time.sleep(np.random.random())
            continue

        jobs        = batch["jobs"]
        params      = batch["params"]
        step        = batch["step"]
        lineage_str = batch["lineage_str"]

        # dvv_lineage = lineage up to but not including the dvv step itself.
        # The parent DTT step is the last entry in that list.
        dvv_lineage   = batch["lineage_names_upstream"]
        root          = params.global_.output_folder
        dvv_step_name = step.step_name

        if not dvv_lineage:
            logger.error(f"Empty lineage for {step_category} job — skipping")
            massive_update_job(db, jobs, "F")
            continue

        # Parent step name is the last element of the upstream lineage
        parent_step_name = dvv_lineage[-1]
        # Parent lineage for reading DTT files = dvv_lineage (same path)
        parent_lineage = dvv_lineage

        logger.info(f"New DVV Job: lineage={lineage_str}")

        # Build (sta1, sta2) pairs as NET.STA.LOC strings,
        # enumerating all configured location code combinations.
        include_single_station = (
            len(params.cc.components_to_compute_single_station) >= 1
            and (params.filter.SC or params.filter.AC)
        )
        all_pairs = []
        for s1, s2 in get_station_pairs(db, include_single_station=include_single_station):
            for loc1 in (s1.locs() or ["00"]):
                for loc2 in (s2.locs() or ["00"]):
                    all_pairs.append((
                        f"{s1.net}.{s1.sta}.{loc1}",
                        f"{s2.net}.{s2.sta}.{loc2}",
                    ))

        # Only include components that were actually computed per filter flags
        cc_comps = [c for c in params.cc.components_to_compute if c]                    if params.filter.CC else []
        ss_comps = [c for c in params.cc.components_to_compute_single_station if c]                    if (params.filter.SC or params.filter.AC) else []
        all_components = list(dict.fromkeys(cc_comps + ss_comps))

        split_pair_type  = params.category_layer.dvv_split_pair_type
        split_components = params.category_layer.dvv_split_components

        # split=Y → one file per individual type/component
        # split=N → one combined ALL file
        if split_pair_type:
            pair_types = []
            if params.filter.CC:
                pair_types.append("CC")
            if params.filter.SC:
                pair_types.append("SC")
            if params.filter.AC:
                pair_types.append("AC")
            if not pair_types:
                pair_types = ["ALL"]  # safety: if all disabled, aggregate everything
        else:
            pair_types = ["ALL"]
        comp_groups = list(all_components) if split_components else ["ALL"]

        mov_stacks = params.stack.mov_stack

        for mov_stack in mov_stacks:
            for pt in pair_types:
                for comp in comp_groups:
                    try:
                        ds_out, ds_pairs = aggregate_dvv_pairs(
                            root=root,
                            parent_lineage=parent_lineage,
                            parent_step_name=parent_step_name,
                            parent_category=parent_category,
                            mov_stack=mov_stack,
                            component=comp,
                            pair_type=pt,
                            pairs=all_pairs,
                            params=params,
                        )
                    except ValueError:
                        logger.debug(
                            f"No data for mov_stack={mov_stack} "
                            f"comp={comp} pair_type={pt} — skipping"
                        )
                        continue

                    try:
                        xr_save_dvv_agg(
                            root, dvv_lineage, dvv_step_name,
                            mov_stack, pt, comp, ds_out,
                        )
                        xr_save_dvv_pairs(
                            root, dvv_lineage, dvv_step_name,
                            mov_stack, pt, comp, ds_pairs,
                        )
                        logger.info(
                            f"Saved dvv_{pt}_{comp} mov_stack={mov_stack} "
                            f"({ds_out.sizes['times']} time steps, "
                            f"{ds_pairs.sizes['pair']} pairs)"
                        )
                    except Exception as e:
                        logger.error(
                            f"Failed saving dvv_{pt}_{comp} "
                            f"mov_stack={mov_stack}: {e}"
                        )

        massive_update_job(db, jobs, "D")

    logger.info(f"*** Finished: {step_category} ***")


if __name__ == "__main__":
    main()
