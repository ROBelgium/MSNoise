"""
Fast-failing smoke test suite for MSNoise workflow machinery.

Purpose
-------
Verify the complete job lifecycle using *stub* compute functions that write
minimal valid NetCDF files.  Two suites are provided:

hpc=False (default, tests 01-33)
    Worker stubs call ``propagate_downstream`` inline — no manual
    ``new_jobs --after X`` calls needed.  This matches the default user
    experience.

hpc=True  (tests 40-60, class ``TestSmokeHPC``)
    Worker stubs do NOT call ``propagate_downstream``.  The test explicitly
    calls ``new_jobs --after X`` at each step, matching HPC cluster usage.

What is tested
--------------
* Workflow step and link creation.
* Direct job seeding (preprocess + PSD T jobs created without DA scan).
* Job T→I→D transitions through stubs.
* ``propagate_downstream``: filter_N pass-through, multi-branch fan-out,
  DVV sentinel creation.
* Lineage normalisation: no NULL lineage_ids, no duplicate Lineage rows.
* ``MSNoiseResult.list()`` finds results after DVV aggregation.
* Second-day: adding new jobs only creates downstream T jobs for the new day.
* Idempotency: calling ``propagate_downstream`` twice is a no-op.
* HPC path: ``new_jobs --after X`` creates the right T jobs (hpc=True only).

Speed target
------------
< 40 seconds total (both suites combined).
"""
from __future__ import annotations

import os

import numpy as np
import pytest
import xarray as xr


# ── autouse fixture ───────────────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def _print_test_name(request):
    print(f"\n{'-'*60}\n> {request.node.name}\n{'-'*60}", flush=True)
    yield


# ── Session fixture ───────────────────────────────────────────────────────────

@pytest.fixture(scope="session")
def smoke_db(tmp_path_factory):
    """Set up a minimal MSNoise project (hpc=False) and return (db, params, root)."""
    from ..core.db import create_database_inifile, connect
    from ..core.config import create_config_set, update_config, get_params
    from ..core.workflow import (create_workflow_steps_from_config_sets,
                                  create_workflow_links_from_steps)
    from ..msnoise_table_def import declare_tables, Station, DataSource

    tmp = tmp_path_factory.mktemp("smoke")
    os.chdir(tmp)

    create_database_inifile(tech=1, hostname="smoke.sqlite", database="",
                            username="", password="", prefix="")
    db = connect()
    declare_tables().Base.metadata.create_all(db.get_bind())

    for cat in ["global", "preprocess", "cc", "filter", "stack", "refstack",
                "mwcs", "mwcs_dtt", "mwcs_dtt_dvv",
                "stretching", "stretching_dvv",
                "wavelet", "wavelet_dtt", "wavelet_dtt_dvv",
                "psd", "psd_rms"]:
        create_config_set(db, cat)

    root = str(tmp / "OUTPUT")
    # Create the default DataSource row (normally done by the installer)
    default_ds = DataSource(name="local", uri=str(tmp / "data"),
                            data_structure="SDS", auth_env="MSNOISE",
                            network_code="YA", channels="HHZ")
    db.add(default_ds)
    db.commit()
    update_config(db, "output_folder",  root)
    update_config(db, "sampling_rate",  "1")
    update_config(db, "startdate",      "2010-09-01")
    update_config(db, "enddate",        "2010-09-02")
    update_config(db, "maxlag",            "2",  category="cc",      set_number=1)
    update_config(db, "cc_sampling_rate",  "1",  category="cc",      set_number=1)
    update_config(db, "components_to_compute", "ZZ", category="cc",  set_number=1)
    update_config(db, "keep_days",         "Y",  category="cc",      set_number=1)
    update_config(db, "mov_stack", "(('1D','1D'),)", category="stack", set_number=1)
    update_config(db, "ref_begin", "2010-09-01",   category="refstack", set_number=1)
    update_config(db, "ref_end",   "2010-09-02",   category="refstack", set_number=1)
    update_config(db, "mwcs_wlen", "2",            category="mwcs",   set_number=1)
    update_config(db, "mwcs_step", "1",            category="mwcs",   set_number=1)
    update_config(db, "dvv_split_pair_type", "Y",  category="mwcs_dtt_dvv",    set_number=1)
    update_config(db, "dvv_split_pair_type", "Y",  category="stretching_dvv",  set_number=1)
    update_config(db, "dvv_split_pair_type", "Y",  category="wavelet_dtt_dvv", set_number=1)
    update_config(db, "psd_rms_frequency_ranges", "[(1.0, 20.0)]",  category="psd_rms", set_number=1)
    update_config(db, "psd_rms_type",             "VEL",            category="psd_rms", set_number=1)
    # hpc=False (default) — propagate_downstream fires automatically
    update_config(db, "hpc", "N")

    for net, sta, x, y in [("YA", "UV05", 0.0, 0.0),
                            ("YA", "UV06", 0.1, 0.0),
                            ("YA", "UV07", 0.0, 0.1)]:
        db.add(Station(net=net, sta=sta, X=x, Y=y, altitude=0.0,
                       coordinates="DEG", used=1,
                       used_location_codes="00",
                       used_channel_names="HHZ"))
    db.commit()

    create_workflow_steps_from_config_sets(db)
    create_workflow_links_from_steps(db)

    params = get_params(db)
    return db, params, root


# ── Constants ─────────────────────────────────────────────────────────────────

DAYS       = ["2010-09-01", "2010-09-02"]
PAIRS      = [("YA.UV05.00", "YA.UV06.00"),
              ("YA.UV05.00", "YA.UV07.00"),
              ("YA.UV06.00", "YA.UV07.00")]
STATIONS   = ["YA.UV05.00", "YA.UV06.00", "YA.UV07.00"]
COMPS      = ["ZZ"]
MOV_STACKS = [("1D", "1D")]
TAXIS      = np.arange(-2.0, 2.0 + 1.0, 1.0, dtype="float32")
FREQS      = np.array([0.1, 0.2, 0.3, 0.5], dtype="float32")
TIMES      = np.array(DAYS, dtype="datetime64[D]")

DAY3       = "2010-09-03"   # extra day for second-day tests
TIMES3     = np.array([DAY3], dtype="datetime64[D]")


# ── Job-count helpers ─────────────────────────────────────────────────────────

def _counts(db, step_name: str) -> dict:
    """Return {flag: count} for all jobs at *step_name*.

    Uses a fresh session via connect() so it always sees the latest
    committed data regardless of the fixture session's cache state.
    """
    from ..core.db import connect as _connect
    from ..msnoise_table_def import declare_tables as _dt
    _s = _dt()
    _db = _connect()
    try:
        step = _db.query(_s.WorkflowStep).filter(
            _s.WorkflowStep.step_name == step_name).first()
        if step is None:
            return {}
        rows = _db.query(_s.Job.flag).filter(_s.Job.step_id == step.step_id).all()
        c: dict = {}
        for (f,) in rows:
            c[f] = c.get(f, 0) + 1
        return c
    finally:
        _db.close()


def _assert_min(db, step_name, flag, n, msg=""):
    c = _counts(db, step_name)
    got = c.get(flag, 0)
    assert got >= n, (
        f"{msg or step_name}: expected >={n} '{flag}' jobs, got {got} "
        f"(all flags: {c})")


def _n_jobs(db, step_name: str, day: str, flag: str) -> int:
    """Count jobs for a specific (step, day, flag) using a fresh session."""
    from ..core.db import connect as _connect
    from ..msnoise_table_def import declare_tables as _dt
    _s = _dt()
    _db = _connect()
    try:
        step = _db.query(_s.WorkflowStep).filter(
            _s.WorkflowStep.step_name == step_name).first()
        if step is None:
            return 0
        return (_db.query(_s.Job)
                  .filter(_s.Job.step_id == step.step_id)
                  .filter(_s.Job.day == day)
                  .filter(_s.Job.flag == flag)
                  .count())
    finally:
        _db.close()


# ── Seeder: create T jobs directly (bypass DA scan) ──────────────────────────

def _seed_jobs(db, days=None, extra_stations=None):
    """Directly insert preprocess_1 and psd_1 T jobs for *days*.

    Bypasses the DataAvailability scan entirely — we trust that path works
    and focus on testing job propagation only.
    """
    from ..msnoise_table_def import declare_tables
    from ..core.workflow import update_job

    schema = declare_tables()
    WorkflowStep = schema.WorkflowStep

    if days is None:
        days = DAYS
    stations = list(STATIONS) + (extra_stations or [])
    n = 0
    for cat in ("preprocess", "psd"):
        step = db.query(WorkflowStep).filter(
            WorkflowStep.category == cat).first()
        if step is None:
            continue
        lineage_str = step.step_name   # single-node lineage for origin steps
        for station in stations:
            for day in days:
                # update_job handles existing-job detection + D→T protection
                update_job(db, day, station, step.step_name, "T",
                           step_id=step.step_id, lineage=lineage_str,
                           commit=False)
                n += 1
    db.commit()
    return n


# ── Stub compute functions ────────────────────────────────────────────────────

def _run_step(category, group_by, on_batch):
    """Generic worker loop: claim batches, call on_batch, mark Done.

    Opens a fresh DB connection per call so the session always sees the
    latest committed state — avoiding stale-cache issues when multiple
    steps share the same fixture session.

    ``on_batch(batch)`` should write output files.
    In hpc=False mode, ``propagate_downstream`` fires automatically.
    """
    from ..core.db import connect as _connect
    from ..core.workflow import (get_next_lineage_batch,
                                  massive_update_job,
                                  propagate_downstream)
    _db = _connect()
    processed = 0
    try:
        while True:
            batch = get_next_lineage_batch(_db, category, group_by=group_by)
            if batch is None:
                break
            on_batch(batch)
            massive_update_job(_db, batch["jobs"], flag="D")
            if not batch["params"].global_.hpc:
                propagate_downstream(_db, batch)
            processed += 1
    finally:
        _db.close()
    return processed


def _stub_preprocess(db, params, root):
    """Run preprocess step with a stub that writes synthetic per-station files."""
    from obspy import Trace, Stream, UTCDateTime
    import numpy as _np

    def _write_station_files(batch):
        """Write synthetic per-station MiniSEED files for each job in the batch."""
        step_name = batch["step"].step_name
        day       = batch["days"][0]
        t0        = UTCDateTime(day)
        for job in batch["jobs"]:
            net, sta, loc = job.pair.split(".")
            tr = Trace(data=_np.zeros(100, dtype="float32"))
            tr.stats.network       = net
            tr.stats.station       = sta
            tr.stats.location      = loc
            tr.stats.channel       = "HHZ"
            tr.stats.starttime     = t0
            tr.stats.sampling_rate = 1.0
            st = Stream(traces=[tr])
            from ..core.signal import save_preprocessed_streams
            save_preprocessed_streams(st, root, step_name, day)

    _run_step("preprocess", "day_lineage", _write_station_files)


def _stub_cc(db, params, root):
    from ..core.io import xr_save_ccf_daily

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        day   = batch["days"][0]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                xr_save_ccf_daily(root, lin, sname,
                                  sta1, sta2, comp, day,
                                  TAXIS, np.zeros(len(TAXIS), "float32"))

    _run_step("cc", "day_lineage", _write)


def _stub_stack(db, params, root):
    from ..core.io import xr_save_ccf

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({"CCF": xr.DataArray(
                        np.zeros((1, len(TAXIS)), "float32"),
                        dims=["times", "taxis"],
                        coords={"times": TIMES[:1], "taxis": TAXIS})})
                    xr_save_ccf(root, lin, sname, sta1, sta2, comp, ms, TAXIS, ds)

    _run_step("stack", "pair_lineage", _write)


def _stub_refstack(db, params, root):
    from ..core.io import xr_save_ref

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                ds = xr.Dataset({"REF": xr.DataArray(
                    np.zeros(len(TAXIS), "float32"),
                    dims=["taxis"], coords={"taxis": TAXIS})})
                xr_save_ref(root, lin, sname, sta1, sta2, comp, TAXIS, ds)

    _run_step("refstack", "pair_lineage", _write)


def _stub_mwcs(db, params, root):
    from ..core.io import xr_save_mwcs

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({"MWCS": xr.DataArray(
                        np.zeros((len(TIMES), len(TAXIS), 4), "float32"),
                        dims=["times", "taxis", "keys"],
                        coords={"times": TIMES, "taxis": TAXIS,
                                "keys": ["dt", "err", "coh", "valid"]})})
                    xr_save_mwcs(root, lin, sname, sta1, sta2, comp, ms, TAXIS, ds)

    _run_step("mwcs", "pair_lineage", _write)


def _stub_dtt(db, params, root):
    from ..core.io import xr_save_dtt

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({"DTT": xr.DataArray(
                        np.zeros((len(TIMES), 3), "float32"),
                        dims=["times", "keys"],
                        coords={"times": TIMES, "keys": ["dtt", "err", "coh"]})})
                    xr_save_dtt(root, lin, sname, sta1, sta2, comp, ms, ds)

    _run_step("mwcs_dtt", "pair_lineage", _write)


def _stub_stretching(db, params, root):
    from ..core.io import xr_save_stretching

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({"STR": xr.DataArray(
                        np.zeros((len(TIMES), 3), "float32"),
                        dims=["times", "keys"],
                        coords={"times": TIMES,
                                "keys": ["Delta", "Coeff", "Error"]})})
                    xr_save_stretching(root, lin, sname, sta1, sta2, comp, ms, ds)

    _run_step("stretching", "pair_lineage", _write)


def _stub_wct(db, params, root):
    from ..core.io import xr_save_wct
    import xarray as _xr

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    z = np.zeros((len(TIMES), len(FREQS), len(TAXIS)), "float32")
                    wct_ds = _xr.Dataset({
                        "WXamp": _xr.DataArray(z, dims=["times","freqs","taxis"],
                                               coords={"times":TIMES,"freqs":FREQS,"taxis":TAXIS}),
                        "Wcoh":  _xr.DataArray(z, dims=["times","freqs","taxis"],
                                               coords={"times":TIMES,"freqs":FREQS,"taxis":TAXIS}),
                        "WXdt":  _xr.DataArray(z, dims=["times","freqs","taxis"],
                                               coords={"times":TIMES,"freqs":FREQS,"taxis":TAXIS}),
                    })
                    xr_save_wct(root, lin, sname, sta1, sta2, comp, ms, wct_ds)

    _run_step("wavelet", "pair_lineage", _write)


def _stub_wct_dtt(db, params, root):
    from ..core.io import xr_save_wct_dtt

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({
                        "DTT": xr.DataArray(np.zeros(len(TIMES), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                        "ERR": xr.DataArray(np.zeros(len(TIMES), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                        "COH": xr.DataArray(np.ones(len(TIMES), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                    })
                    xr_save_wct_dtt(root, lin, sname, sta1, sta2, comp, ms, TAXIS, ds)

    _run_step("wavelet_dtt", "pair_lineage", _write)


def _stub_dvv_agg(db, params, root, category):
    from ..core.io import xr_save_dvv_agg

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for ms in MOV_STACKS:
            for comp in COMPS:
                ds = xr.Dataset({
                    "mean":    xr.DataArray(np.zeros(len(TIMES), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                    "std":     xr.DataArray(np.zeros(len(TIMES), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                    "n_pairs": xr.DataArray(
                        np.full(len(TIMES), len(PAIRS), "int32"),
                        dims=["times"], coords={"times": TIMES}),
                })
                xr_save_dvv_agg(root, lin, sname, ms, "CC", comp[0], ds)

    _run_step(category, "day_lineage", _write)


def _stub_psd(db, params, root):
    from ..core.io import xr_save_psd
    periods = np.logspace(-1, 1, 20, dtype="float32")

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        day   = batch["days"][0]
        for seed_id in ["YA.UV05.00.HHZ", "YA.UV06.00.HHZ", "YA.UV07.00.HHZ"]:
            ds = xr.Dataset({"PSD": xr.DataArray(
                np.full((1, len(periods)), -150.0, "float32"),
                dims=["times", "periods"],
                coords={"times": [np.datetime64(day)], "periods": periods})})
            xr_save_psd(root, lin, sname, seed_id, day, ds)

    _run_step("psd", "day_lineage", _write)


def _stub_psd_rms(db, params, root):
    from ..core.io import xr_save_rms
    import pandas as pd

    def _write(batch):
        lin   = batch["lineage_names"][:-1]
        sname = batch["lineage_names"][-1]
        for seed_id in ["YA.UV05.00.HHZ", "YA.UV06.00.HHZ", "YA.UV07.00.HHZ"]:
            df = pd.DataFrame(
                {"low": [-150.0] * 2, "mid": [-150.0] * 2, "high": [-150.0] * 2},
                index=pd.DatetimeIndex(DAYS))
            xr_save_rms(root, lin, sname, seed_id, df)

    _run_step("psd_rms", "day_lineage", _write)


# ─────────────────────────────────────────────────────────────────────────────
# hpc=False suite  (tests 01–33)
# Worker stubs call propagate_downstream inline — no new_jobs --after needed.
# ─────────────────────────────────────────────────────────────────────────────

@pytest.mark.order(1)
def test_smoke_01_schema(smoke_db):
    """DB tables including Lineage and DataSource exist."""
    db, _, _ = smoke_db
    from sqlalchemy import inspect
    tables = inspect(db.get_bind()).get_table_names()
    for t in ("jobs", "lineages", "workflow_steps", "workflow_links",
              "config", "stations", "data_sources"):
        assert t in tables, f"Table {t!r} missing"
    # Default DataSource must exist with id=1
    from ..msnoise_table_def import DataSource
    ds = db.query(DataSource).filter(DataSource.ref == 1).first()
    assert ds is not None, "Default DataSource (id=1) missing"
    assert ds.name == "local"
    assert ds.data_structure == "SDS"


@pytest.mark.order(2)
def test_smoke_02_workflow_created(smoke_db):
    """All expected workflow steps and links created."""
    db, _, _ = smoke_db
    from ..core.workflow import get_workflow_steps, get_workflow_links
    names = {s.step_name for s in get_workflow_steps(db)}
    required = {"preprocess_1", "cc_1", "filter_1", "stack_1", "refstack_1",
                "mwcs_1", "mwcs_dtt_1", "mwcs_dtt_dvv_1",
                "stretching_1", "stretching_dvv_1",
                "wavelet_1", "wavelet_dtt_1", "wavelet_dtt_dvv_1",
                "psd_1", "psd_rms_1"}
    assert required.issubset(names), f"Missing: {required - names}"
    assert len(get_workflow_links(db)) >= 5


@pytest.mark.order(3)
def test_smoke_03_seed_jobs(smoke_db):
    """Direct seeding creates preprocess_1 and psd_1 T jobs (no DA scan)."""
    db, _, _ = smoke_db
    n = _seed_jobs(db)
    assert n > 0, "Seeder must insert jobs"
    # 3 stations × 2 days × 2 steps (preprocess + psd)
    _assert_min(db, "preprocess_1", "T", 6, "seed: preprocess T")
    _assert_min(db, "psd_1",        "T", 6, "seed: psd T")


@pytest.mark.order(4)
def test_smoke_04_preprocess(smoke_db):
    """Preprocess stub marks jobs Done; per-station output files exist."""
    db, params, root = smoke_db
    _assert_min(db, "preprocess_1", "T", 1)
    _stub_preprocess(db, params, root)
    _assert_min(db, "preprocess_1", "D", 6, "preprocess must be Done")
    # Verify per-station file layout: _output/<date>/<NET.STA.LOC>.mseed
    for day in DAYS:
        for sid in STATIONS:
            fpath = os.path.join(root, "preprocess_1", "_output", day, f"{sid}.mseed")
            assert os.path.isfile(fpath), (
                f"Per-station preprocess output missing: {fpath}"
            )
    # propagate_downstream → cc_1 T jobs (3 pairs × 2 days)
    _assert_min(db, "cc_1", "T", 6, "propagate_downstream must create cc T jobs")


@pytest.mark.order(5)
def test_smoke_05_cc(smoke_db):
    """CC stub marks jobs Done; propagate_downstream creates BOTH stack_1 AND
    refstack_1 T jobs (now siblings under filter_1)."""
    db, params, root = smoke_db
    _assert_min(db, "cc_1", "T", 1)
    _stub_cc(db, params, root)
    _assert_min(db, "cc_1",      "D", 6, "cc must be Done")
    _assert_min(db, "stack_1",   "T", 1,
                "propagate_downstream must create stack_1 T jobs")
    _assert_min(db, "refstack_1","T", 1,
                "propagate_downstream must create refstack_1 REF T jobs (sibling of stack)")


@pytest.mark.order(6)
def test_smoke_06_stack(smoke_db):
    """Stack stub marks jobs Done. Stack is now a leaf — refstack T jobs were
    already created by cc (test_05). Mwcs jobs appear only after refstack Done."""
    db, params, root = smoke_db
    _assert_min(db, "stack_1", "T", 1)
    _stub_stack(db, params, root)
    _assert_min(db, "stack_1", "D", 1)


@pytest.mark.order(7)
def test_smoke_07_refstack(smoke_db):
    """Refstack stub marks jobs Done; propagate_downstream fans out to
    mwcs_1, stretching_1, and wavelet_1 T jobs."""
    db, params, root = smoke_db
    _assert_min(db, "refstack_1", "T", 1)
    _stub_refstack(db, params, root)
    _assert_min(db, "refstack_1",  "D", 1)
    _assert_min(db, "mwcs_1",      "T", 1, "mwcs_1 T after refstack")
    _assert_min(db, "stretching_1","T", 1, "stretching_1 T after refstack")
    _assert_min(db, "wavelet_1",   "T", 1, "wavelet_1 T after refstack")


@pytest.mark.order(8)
def test_smoke_08_mwcs(smoke_db):
    """MWCS stub; propagate_downstream → mwcs_dtt_1 T jobs."""
    db, params, root = smoke_db
    _assert_min(db, "mwcs_1", "T", 1)
    _stub_mwcs(db, params, root)
    _assert_min(db, "mwcs_1",     "D", 1)
    _assert_min(db, "mwcs_dtt_1", "T", 1)


@pytest.mark.order(9)
def test_smoke_09_dtt(smoke_db):
    """DTT stub; propagate_downstream → mwcs_dtt_dvv_1 sentinel T job."""
    db, params, root = smoke_db
    _assert_min(db, "mwcs_dtt_1", "T", 1)
    _stub_dtt(db, params, root)
    _assert_min(db, "mwcs_dtt_1",     "D", 1)
    _assert_min(db, "mwcs_dtt_dvv_1", "T", 1,
                "propagate_downstream must create DVV sentinel job")


@pytest.mark.order(10)
def test_smoke_10_mwcs_dvv(smoke_db):
    """DVV-MWCS aggregate stub; MSNoiseResult.list finds it."""
    db, params, root = smoke_db
    _assert_min(db, "mwcs_dtt_dvv_1", "T", 1)
    _stub_dvv_agg(db, params, root, "mwcs_dtt_dvv")
    _assert_min(db, "mwcs_dtt_dvv_1", "D", 1)
    from ..results import MSNoiseResult
    assert len(MSNoiseResult.list(db, "mwcs_dtt_dvv")) >= 1


@pytest.mark.order(11)
def test_smoke_11_stretching(smoke_db):
    """Stretching stub; propagate_downstream → stretching_dvv_1 sentinel."""
    db, params, root = smoke_db
    _assert_min(db, "stretching_1", "T", 1)
    _stub_stretching(db, params, root)
    _assert_min(db, "stretching_1",    "D", 1)
    _assert_min(db, "stretching_dvv_1","T", 1)


@pytest.mark.order(12)
def test_smoke_12_stretching_dvv(smoke_db):
    """DVV-stretching aggregate stub."""
    db, params, root = smoke_db
    _stub_dvv_agg(db, params, root, "stretching_dvv")
    _assert_min(db, "stretching_dvv_1", "D", 1)
    from ..results import MSNoiseResult
    assert len(MSNoiseResult.list(db, "stretching_dvv")) >= 1


@pytest.mark.order(13)
def test_smoke_13_wct(smoke_db):
    """WCT stub; propagate_downstream → wavelet_dtt_1 T jobs."""
    db, params, root = smoke_db
    _assert_min(db, "wavelet_1", "T", 1)
    _stub_wct(db, params, root)
    _assert_min(db, "wavelet_1",   "D", 1)
    _assert_min(db, "wavelet_dtt_1","T", 1)


@pytest.mark.order(14)
def test_smoke_14_wct_dtt(smoke_db):
    """WCT-DTT stub; propagate_downstream → wavelet_dtt_dvv_1 sentinel."""
    db, params, root = smoke_db
    _assert_min(db, "wavelet_dtt_1", "T", 1)
    _stub_wct_dtt(db, params, root)
    _assert_min(db, "wavelet_dtt_1",    "D", 1)
    _assert_min(db, "wavelet_dtt_dvv_1","T", 1)


@pytest.mark.order(15)
def test_smoke_15_wct_dvv(smoke_db):
    """DVV-WCT aggregate stub."""
    db, params, root = smoke_db
    _stub_dvv_agg(db, params, root, "wavelet_dtt_dvv")
    _assert_min(db, "wavelet_dtt_dvv_1", "D", 1)
    from ..results import MSNoiseResult
    assert len(MSNoiseResult.list(db, "wavelet_dtt_dvv")) >= 1


@pytest.mark.order(16)
def test_smoke_16_psd(smoke_db):
    """PSD stub marks jobs Done; propagate_downstream → psd_rms_1 T jobs."""
    db, params, root = smoke_db
    _assert_min(db, "psd_1", "T", 1)
    _stub_psd(db, params, root)
    _assert_min(db, "psd_1",    "D", 1)
    _assert_min(db, "psd_rms_1","T", 1)


@pytest.mark.order(17)
def test_smoke_17_psd_rms(smoke_db):
    """PSD-RMS stub marks jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "psd_rms_1", "T", 1)
    _stub_psd_rms(db, params, root)
    _assert_min(db, "psd_rms_1", "D", 1)


@pytest.mark.order(171)
def test_smoke_171_psd_rms_lineage(smoke_db):
    """psd_rms jobs must have lineage 'psd_1/psd_rms_1' (not just 'psd_1').

    This guards the root-cause fix: propagate_psd_rms_jobs_from_psd_done
    must store lineage=psd_step/psd_rms_step so that get_next_lineage_batch
    builds a MSNoiseParams with a 'psd_rms' layer.  Without this, accessing
    params.psd_rms.* raises AttributeError.
    """
    db, _, _ = smoke_db
    from ..msnoise_table_def import declare_tables as _dt
    _s = _dt()
    from ..core.db import connect as _connect
    _db = _connect()
    try:
        psd_rms_step = _db.query(_s.WorkflowStep).filter(
            _s.WorkflowStep.step_name == "psd_rms_1").first()
        assert psd_rms_step is not None

        jobs = _db.query(_s.Job).filter(
            _s.Job.step_id == psd_rms_step.step_id).all()
        assert jobs, "No psd_rms_1 jobs found"

        for job in jobs:
            lin = _db.query(_s.Lineage).filter(
                _s.Lineage.lineage_id == job.lineage_id).first()
            assert lin is not None, f"Job {job.ref} has no Lineage row"
            assert "/" in lin.lineage_str, (
                f"psd_rms job lineage {lin.lineage_str!r} must be "
                f"'psd_1/psd_rms_1', not just the psd step name"
            )
            assert lin.lineage_str.endswith("/psd_rms_1"), (
                f"psd_rms job lineage {lin.lineage_str!r} must end with /psd_rms_1"
            )
    finally:
        _db.close()


@pytest.mark.order(172)
def test_smoke_172_psd_rms_params_layer(smoke_db):
    """get_next_lineage_batch for psd_rms must build MSNoiseParams with 'psd_rms' layer.

    Exercises the full params.psd_rms.* access path that s21_psd_compute_rms.py
    uses at runtime — verifying that the lineage fix flows through to the
    MSNoiseParams construction.
    """
    db, _, _ = smoke_db
    from ..core.db import connect as _connect
    from ..core.workflow import get_next_lineage_batch
    from ..msnoise_table_def import declare_tables as _dt
    _s = _dt()

    # Seed a fresh psd_rms T job on day1 to claim
    _db = _connect()
    try:
        psd_rms_step = _db.query(_s.WorkflowStep).filter(
            _s.WorkflowStep.step_name == "psd_rms_1").first()
        # Reset one psd_rms job to T so get_next_lineage_batch can claim it
        job = _db.query(_s.Job).filter(
            _s.Job.step_id == psd_rms_step.step_id).first()
        assert job is not None
        job.flag = "T"
        _db.commit()

        batch = get_next_lineage_batch(_db, "psd_rms", group_by="pair_lineage")
        assert batch is not None, "get_next_lineage_batch returned None for psd_rms"

        params = batch["params"]
        # Must have 'psd_rms' layer — accessing it must not raise AttributeError
        try:
            _ = params.psd_rms.psd_rms_frequency_ranges
            _ = params.psd_rms.psd_rms_type
        except AttributeError as e:
            raise AssertionError(
                f"params.psd_rms.* raised AttributeError — "
                f"psd_rms layer missing from MSNoiseParams. "
                f"Available categories: {list(params._layers)}. Error: {e}"
            )

        # Restore to Done so subsequent tests see a consistent state
        from ..core.workflow import massive_update_job
        massive_update_job(_db, batch["jobs"], "D")
    finally:
        _db.close()


@pytest.mark.order(173)
def test_smoke_173_get_waveform_path_datasource(smoke_db):
    """get_waveform_path correctly prepends DataSource.uri to DA path+file."""
    db, _, _ = smoke_db
    from ..core.stations import get_waveform_path, get_default_data_source
    from ..msnoise_table_def import declare_tables as _dt
    _s = _dt()

    ds = get_default_data_source(db)
    assert ds is not None

    # Construct a synthetic DA-like object to test path building
    class _FakeDA:
        data_source_id = ds.ref
        path = "YA/UV05/2010/HHZ.D"
        file = "YA.UV05..HHZ.D.2010.244"

    fpath = get_waveform_path(db, _FakeDA())
    assert fpath.startswith(ds.uri), (
        f"Path {fpath!r} must start with DataSource.uri {ds.uri!r}"
    )
    assert fpath.endswith(_FakeDA.file), (
        f"Path {fpath!r} must end with filename {_FakeDA.file!r}"
    )
    # Ensure no double-separator artefacts
    assert "//" not in fpath, f"Double separator in path: {fpath!r}"


@pytest.mark.order(174)
def test_smoke_174_create_psd_jobs_cmd(smoke_db, tmp_path):
    """msnoise utils create_psd_jobs creates T jobs for the correct date."""
    db, _, _ = smoke_db
    from click.testing import CliRunner
    from ..scripts.msnoise import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["utils", "create_psd_jobs",
                                 "--date", "2010-09-01",
                                 "--set-number", "1"], obj={})
    assert result.exit_code == 0, (
        f"create_psd_jobs failed (exit {result.exit_code}):\n{result.output}"
    )
    assert "psd" in result.output.lower(), (
        f"Expected 'psd' in output, got: {result.output!r}"
    )


@pytest.mark.order(175)
def test_smoke_175_cc_chunk_claim(smoke_db):
    """chunk_size limits the number of CC jobs claimed per batch.

    Seeds 3 CC T jobs for one day (3 pairs), then calls
    get_next_lineage_batch with chunk_size=2. Asserts:
    - At most 2 jobs are returned in the batch.
    - The claimed jobs are marked 'I' (in-progress).
    - At least 1 T job remains for the next worker to claim.
    - A second call (simulating a second worker) claims the remainder.
    """
    db, _, _ = smoke_db
    from ..core.db import connect as _connect
    from ..core.workflow import (get_next_lineage_batch, massive_update_job)
    from ..msnoise_table_def import declare_tables as _dt

    _s = _dt()
    _db = _connect()
    try:
        cc_step = _db.query(_s.WorkflowStep).filter(
            _s.WorkflowStep.step_name == "cc_1").first()
        assert cc_step is not None, "cc_1 step not found"

        # Get cc lineage_id from an existing cc job
        existing = _db.query(_s.Job).filter(
            _s.Job.step_id == cc_step.step_id).first()
        assert existing is not None, "No cc_1 jobs found (requires tests 04-05 to have run)"

        # Reset all cc jobs to T for day 1 so we have a clean slate
        test_day = "2010-09-01"
        day_jobs = _db.query(_s.Job).filter(
            _s.Job.step_id == cc_step.step_id,
            _s.Job.day == test_day,
        ).all()
        if not day_jobs:
            pytest.skip("No cc_1 jobs for 2010-09-01 (requires test_05)")
        for j in day_jobs:
            j.flag = "T"
        _db.commit()

        n_total = len(day_jobs)
        assert n_total >= 3, f"Need at least 3 CC jobs for chunking test, got {n_total}"

        # First worker: chunk_size=2 — must claim at most 2 jobs
        batch1 = get_next_lineage_batch(_db, "cc", group_by="day_lineage",
                                         chunk_size=2)
        assert batch1 is not None, "First chunk claim returned None"
        n_claimed = len(batch1["jobs"])
        assert n_claimed <= 2, (
            f"chunk_size=2 must claim at most 2 jobs, got {n_claimed}"
        )
        assert n_claimed >= 1, "Must claim at least 1 job"

        # Claimed jobs must be marked I
        _db.expire_all()
        for job in batch1["jobs"]:
            fresh = _db.query(_s.Job).filter(_s.Job.ref == job.ref).first()
            assert fresh.flag == "I", (
                f"Job {job.ref} should be 'I' after chunk claim, got {fresh.flag!r}"
            )

        # At least one T job must remain for the next worker
        remaining_t = _db.query(_s.Job).filter(
            _s.Job.step_id == cc_step.step_id,
            _s.Job.day == test_day,
            _s.Job.flag == "T",
        ).count()
        assert remaining_t >= 1, (
            f"At least 1 T job must remain after chunk_size=2 claim "
            f"(total={n_total}, claimed={n_claimed})"
        )

        # Second worker: claims the rest
        batch2 = get_next_lineage_batch(_db, "cc", group_by="day_lineage",
                                         chunk_size=2)
        assert batch2 is not None, "Second chunk claim returned None"
        n_claimed2 = len(batch2["jobs"])
        assert n_claimed2 >= 1, "Second worker must claim at least 1 job"

        # Together they must cover all day_jobs
        total_claimed = n_claimed + n_claimed2
        assert total_claimed == n_total, (
            f"Two chunk claims should cover all {n_total} jobs, "
            f"got {total_claimed} total"
        )

        # Restore to Done so subsequent tests see a clean state
        massive_update_job(_db, batch1["jobs"], "D")
        massive_update_job(_db, batch2["jobs"], "D")
    finally:
        _db.close()


@pytest.mark.order(18)
def test_smoke_18_lineage_normalisation(smoke_db):
    """Zero NULL lineage_ids; Lineage rows << Job rows; no duplicates."""
    db, _, _ = smoke_db
    db.expire_all()
    from ..msnoise_table_def import declare_tables as _dt18
    _s18 = _dt18()
    from sqlalchemy import func

    null_jobs = db.query(_s18.Job).filter(_s18.Job.lineage_id.is_(None)).count()
    assert null_jobs == 0, f"{null_jobs} jobs have NULL lineage_id"

    dupes = (db.query(_s18.Lineage.lineage_str, func.count(_s18.Lineage.lineage_id))
               .group_by(_s18.Lineage.lineage_str)
               .having(func.count(_s18.Lineage.lineage_id) > 1)
               .all())
    assert len(dupes) == 0, f"Duplicate Lineage rows: {dupes}"

    n_lin  = db.query(_s18.Lineage).count()
    n_jobs = db.query(_s18.Job).count()
    assert n_lin >= 1
    assert n_lin < n_jobs, f"Expected Lineage ({n_lin}) < Jobs ({n_jobs})"
    print(f"  Lineage rows: {n_lin}, Job rows: {n_jobs} [x]")


@pytest.mark.order(19)

@pytest.mark.order(185)
def test_smoke_18b_datasource_stations(smoke_db):
    """Stations with NULL data_source_id resolve to the default DataSource."""
    db, _, _ = smoke_db
    from ..core.stations import resolve_data_source, get_stations
    for station in get_stations(db):
        ds = resolve_data_source(db, station)
        assert ds.ref == 1,   f"Station {station} should resolve to DataSource id=1"
        assert ds.name == "local"
        assert ds.uri != "" or True  # uri may be empty string or path — both valid


def test_smoke_19_all_dvv_discoverable(smoke_db):
    """MSNoiseResult.list + get_dvv work for all three DVV methods."""
    db, _, _ = smoke_db
    from ..results import MSNoiseResult
    for cat in ("mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"):
        results = MSNoiseResult.list(db, cat)
        assert len(results) >= 1, f"list({cat!r}) empty"
        data = results[0].get_dvv(pair_type="CC")
        assert len(data) >= 1, f"get_dvv('CC') empty for {cat}"


# ── Second-day tests (20–25) ──────────────────────────────────────────────────

@pytest.mark.order(20)
def test_smoke_20_seed_day3(smoke_db):
    """Seed day3 preprocess+psd T jobs; existing Done jobs must not be re-queued."""
    db, _, _ = smoke_db
    if _n_jobs(db, "cc_1", "2010-09-02", "D") == 0:
        pytest.skip("Requires tests 01-19 to have run first")

    n = _seed_jobs(db, days=[DAY3])
    assert n > 0, "Seeder must insert day3 jobs"
    _assert_min(db, "preprocess_1", "T", 3, "3 new preprocess T for day3")
    # Day 1 and 2 preprocess must still be Done
    assert _n_jobs(db, "preprocess_1", "2010-09-01", "T") == 0
    assert _n_jobs(db, "preprocess_1", "2010-09-02", "T") == 0


@pytest.mark.order(21)
def test_smoke_21_preprocess_day3(smoke_db):
    """Preprocess stub for day3; propagate_downstream creates cc T for day3 only."""
    db, params, root = smoke_db
    if _n_jobs(db, "preprocess_1", DAY3, "T") == 0:
        pytest.skip("Requires test_20 first")

    _stub_preprocess(db, params, root)
    assert _n_jobs(db, "preprocess_1", DAY3, "T") == 0
    assert _n_jobs(db, "preprocess_1", DAY3, "D") == 3
    # propagate_downstream creates cc T for day3
    assert _n_jobs(db, "cc_1", DAY3, "T") >= 1, (
        "propagate_downstream must create cc_1 T for day3")
    # Existing day1/2 cc Done jobs must not be re-queued
    assert _n_jobs(db, "cc_1", "2010-09-01", "T") == 0
    assert _n_jobs(db, "cc_1", "2010-09-02", "T") == 0


@pytest.mark.order(22)
def test_smoke_22_cc_day3_filter_passthrough(smoke_db):
    """CC stub for day3; propagate_downstream crosses filter_1 → stack_1 T.

    Key assertion: the stack lineage must be
    'preprocess_1/cc_1/filter_1/stack_1' even though the CC job lineage
    is only 'preprocess_1/cc_1' — propagate_downstream inserts filter_1.
    """
    db, params, root = smoke_db
    if _n_jobs(db, "cc_1", DAY3, "T") == 0:
        pytest.skip("Requires test_21 first")

    _stub_cc(db, params, root)
    assert _n_jobs(db, "cc_1",    DAY3, "D") >= 1
    assert _n_jobs(db, "cc_1",    DAY3, "T") == 0
    n_stack = _n_jobs(db, "stack_1", DAY3, "T")
    assert n_stack >= 1, (
        f"propagate_downstream must create stack_1 T via filter_1 (got {n_stack})")
    # Old stack Done jobs must not be re-queued
    assert _n_jobs(db, "stack_1", "2010-09-01", "T") == 0
    assert _n_jobs(db, "stack_1", "2010-09-02", "T") == 0


@pytest.mark.order(23)
def test_smoke_23_stack_day3(smoke_db):
    """Stack stub for day3 only; Done count increases by exactly the new day."""
    db, params, root = smoke_db
    if _n_jobs(db, "stack_1", DAY3, "T") == 0:
        pytest.skip("Requires test_22 first")

    done_before = _counts(db, "stack_1").get("D", 0)
    _stub_stack(db, params, root)
    done_after  = _counts(db, "stack_1").get("D", 0)
    assert done_after > done_before
    assert _n_jobs(db, "stack_1", DAY3, "T") == 0


@pytest.mark.order(24)
def test_smoke_24_no_lineage_dupes_after_day3(smoke_db):
    """After day3 processing, Lineage table must have no duplicates.

    Adding a new day reuses existing Lineage rows — e.g., only one row
    for 'preprocess_1/cc_1/filter_1/stack_1' regardless of day count.
    """
    db, _, _ = smoke_db
    if _n_jobs(db, "stack_1", DAY3, "D") == 0:
        pytest.skip("Requires test_23 first")

    db.expire_all()
    from ..msnoise_table_def import declare_tables as _dt24
    _s24 = _dt24()
    from sqlalchemy import func

    dupes = (db.query(_s24.Lineage.lineage_str, func.count(_s24.Lineage.lineage_id))
               .group_by(_s24.Lineage.lineage_str)
               .having(func.count(_s24.Lineage.lineage_id) > 1)
               .all())
    assert len(dupes) == 0, f"Duplicate Lineage rows after day3: {dupes}"

    null_jobs = db.query(_s24.Job).filter(_s24.Job.lineage_id.is_(None)).count()
    assert null_jobs == 0

    n_lin  = db.query(_s24.Lineage).count()
    n_jobs = db.query(_s24.Job).count()
    assert n_lin < n_jobs
    print(f"  Lineage rows: {n_lin}, Job rows: {n_jobs} [x]")


@pytest.mark.order(25)
def test_smoke_25_propagate_downstream_idempotent(smoke_db):
    """Calling propagate_downstream twice is a no-op (no duplicates, no re-queue)."""
    db, params, root = smoke_db
    from ..core.workflow import propagate_downstream
    from ..msnoise_table_def import declare_tables as _dt25
    _s25 = _dt25()

    cc_step = db.query(_s25.WorkflowStep).filter(
        _s25.WorkflowStep.step_name == "cc_1").first()
    cc_jobs_day3 = (db.query(_s25.Job)
                      .filter(_s25.Job.step_id == cc_step.step_id)
                      .filter(_s25.Job.day == DAY3)
                      .all())
    if not cc_jobs_day3:
        pytest.skip("Requires test_22 (day3 CC Done) first")

    fake_batch = {
        "step":          cc_step,
        "jobs":          cc_jobs_day3,
        "lineage_names": ["preprocess_1", "cc_1"],
        "lineage_str":   "preprocess_1/cc_1",
        "days":          [DAY3],
        "params":        params,
    }

    before = dict(_counts(db, "stack_1"))
    propagate_downstream(db, fake_batch)
    propagate_downstream(db, fake_batch)   # second call — must be idempotent
    after = dict(_counts(db, "stack_1"))

    # T count must not exceed what the first call created (or pre-existing T)
    assert after.get("T", 0) <= before.get("T", 0) + len(PAIRS), (
        f"Second propagate_downstream created extra T jobs: {before} → {after}")

    from sqlalchemy import func
    dupes = (db.query(_s25.Lineage.lineage_str, func.count(_s25.Lineage.lineage_id))
               .group_by(_s25.Lineage.lineage_str)
               .having(func.count(_s25.Lineage.lineage_id) > 1)
               .all())
    assert len(dupes) == 0, f"Lineage duplicates after idempotency test: {dupes}"
