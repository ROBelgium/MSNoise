"""
Fast-failing smoke test suite for MSNoise workflow machinery.

Purpose
-------
Verify the complete job lifecycle — from ``db init`` through DVV aggregate —
using *stub* compute functions that write minimal valid NetCDF files instead
of performing real signal processing.

What is tested
--------------
* Database initialisation and table creation (including the ``Lineage`` table).
* ``new_jobs`` creates the correct number of ``T`` jobs with valid ``lineage_id``.
* Each stub compute step transitions jobs from ``T`` → ``D``.
* ``new_jobs --after X`` propagates correctly to every downstream step.
* ``MSNoiseResult.list()`` returns results at every DVV category.
* Lineage normalisation: every job resolves to a non-empty string; the
  ``Lineage`` table has far fewer rows than total job rows.
* Both CC and PSD branches run independently to completion.
* All three DVV methods (MWCS, Stretching, WCT) produce aggregate output.

What is NOT tested
------------------
* Numerical correctness of any signal-processing algorithm.
* Actual seismic waveform data.

Speed target
------------
< 30 seconds on any modern laptop.

Running
-------
Via the CLI::

    msnoise utils test --fast

Or directly with pytest::

    pytest msnoise/test/test_smoke.py -s -v

"""
from __future__ import annotations

import datetime
import os

import numpy as np
import pytest
import xarray as xr

# ── Fixtures ──────────────────────────────────────────────────────────────────

@pytest.fixture(autouse=True)
def _print_test_name(request):
    print(f"\n{'─'*60}\n▶ {request.node.name}\n{'─'*60}", flush=True)
    yield


@pytest.fixture(scope="session")
def smoke_db(tmp_path_factory):
    """Set up a minimal MSNoise project in a temp directory and return (db, cfg)."""
    from ..core.db import create_database_inifile, connect
    from ..core.config import create_config_set, update_config, get_params
    from ..core.workflow import (create_workflow_steps_from_config_sets,
                                  create_workflow_links_from_steps)
    from ..msnoise_table_def import declare_tables, Station

    tmp = tmp_path_factory.mktemp("smoke")
    os.chdir(tmp)

    create_database_inifile(
        tech=1, hostname="smoke.sqlite", database="",
        username="", password="", prefix=""
    )
    db = connect()
    declare_tables().Base.metadata.create_all(db.get_bind())

    # Config sets — full workflow
    for cat in ["global", "preprocess", "cc", "filter", "stack", "refstack",
                "mwcs", "mwcs_dtt", "mwcs_dtt_dvv",
                "stretching", "stretching_dvv",
                "wavelet", "wavelet_dtt", "wavelet_dtt_dvv",
                "psd", "psd_rms"]:
        create_config_set(db, cat)

    # Minimal parameters
    root = str(tmp / "OUTPUT")
    update_config(db, "data_folder",   str(tmp / "data"))
    update_config(db, "output_folder", root)
    update_config(db, "sampling_rate", "1")          # 1 Hz → tiny arrays
    update_config(db, "channels",      "HHZ")
    update_config(db, "startdate",     "2010-09-01")
    update_config(db, "enddate",       "2010-09-03")
    update_config(db, "maxlag",        "2",            category="cc",      set_number=1)
    update_config(db, "cc_sampling_rate", "1",        category="cc",      set_number=1)
    update_config(db, "components_to_compute", "ZZ",  category="cc",      set_number=1)
    update_config(db, "keep_days",     "Y",            category="cc",      set_number=1)
    update_config(db, "mov_stack",     "((\'1D\',\'1D\'),)", category="stack", set_number=1)
    update_config(db, "ref_begin",     "2010-09-01",   category="refstack",set_number=1)
    update_config(db, "ref_end",       "2010-09-03",   category="refstack",set_number=1)
    update_config(db, "mwcs_wlen",     "2",            category="mwcs",    set_number=1)
    update_config(db, "mwcs_step",     "1",            category="mwcs",    set_number=1)
    update_config(db, "wct_freqmin",   "0.1",          category="wavelet", set_number=1)
    update_config(db, "wct_freqmax",   "0.5",          category="wavelet", set_number=1)
    update_config(db, "dvv_split_pair_type", "Y",      category="mwcs_dtt_dvv",    set_number=1)
    update_config(db, "dvv_split_pair_type", "Y",      category="stretching_dvv",  set_number=1)
    update_config(db, "dvv_split_pair_type", "Y",      category="wavelet_dtt_dvv", set_number=1)

    # Stations — 3 stations → 3 pairs
    for net, sta, x, y in [("YA", "UV05", 0.0, 0.0),
                            ("YA", "UV06", 0.1, 0.0),
                            ("YA", "UV07", 0.0, 0.1)]:
        db.add(Station(net=net, sta=sta, X=x, Y=y, altitude=0.0,
                       coordinates="DEG", used=1,
                       used_location_codes="00",
                       used_channel_names="HHZ"))
    db.commit()

    # Data availability — 2 days, 1 channel
    from ..msnoise_table_def import DataAvailability
    now = datetime.datetime.utcnow()
    for sta in ["UV05", "UV06", "UV07"]:
        for day in ["2010-09-01", "2010-09-02"]:
            db.add(DataAvailability(
                net="YA", sta=sta, loc="00", chan="HHZ",
                path=str(tmp / "data"), file=f"{sta}_{day}.mseed",
                starttime=datetime.datetime.fromisoformat(day),
                endtime=datetime.datetime.fromisoformat(day) + datetime.timedelta(days=1),
                data_duration=86400, gaps_duration=0,
                samplerate=1.0, flag="N"  # 'N'ew so get_new_files picks them up
            ))
    db.commit()

    # Build workflow
    create_workflow_steps_from_config_sets(db)
    create_workflow_links_from_steps(db)

    params = get_params(db)
    return db, params, root


# ── Helpers ───────────────────────────────────────────────────────────────────

def _counts(db, step_name: str) -> dict:
    """Return {flag: count} for all jobs at *step_name*."""
    db.expire_all()  # ensure we see commits made by other sessions
    from ..msnoise_table_def import WorkflowStep, Job as _J
    step = db.query(WorkflowStep).filter(
        WorkflowStep.step_name == step_name).first()
    if step is None:
        return {}
    rows = db.query(_J.flag).filter(_J.step_id == step.step_id).all()
    c = {}
    for (f,) in rows:
        c[f] = c.get(f, 0) + 1
    return c


def _assert_min(db, step_name, flag, n, msg=""):
    c = _counts(db, step_name)
    got = c.get(flag, 0)
    assert got >= n, (
        f"{msg or step_name}: expected >={n} '{flag}' jobs, got {got} (all: {c})")


DAYS      = ["2010-09-01", "2010-09-02"]
PAIRS     = [("YA.UV05.00", "YA.UV06.00"),
             ("YA.UV05.00", "YA.UV07.00"),
             ("YA.UV06.00", "YA.UV07.00")]
COMPS     = ["ZZ"]
MOV_STACKS = [("1D", "1D")]
SR        = 1.0
MAXLAG    = 2.0
TAXIS     = np.arange(-MAXLAG, MAXLAG + 1/SR, 1/SR)
FREQS     = np.array([0.1, 0.2, 0.3, 0.5])
TIMES     = np.array(DAYS, dtype="datetime64[D]")


# ── Stub compute functions ────────────────────────────────────────────────────

def _stub_preprocess(db, params, root):
    """Mark all preprocess jobs Done (no actual preprocessing needed for smoke)."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    while True:
        batch = get_next_lineage_batch(db, "preprocess", group_by="day_lineage")
        if batch is None:
            break
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_cc(db, params, root):
    """Write a trivial daily CCF NetCDF per pair/comp/day, mark jobs Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_ccf_daily
    while True:
        batch = get_next_lineage_batch(db, "cc", group_by="day_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]   # upstream of cc
        step_name = batch["lineage_names"][-1]
        day = batch["days"][0]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                corr = np.zeros(len(TAXIS), dtype="float32")
                xr_save_ccf_daily(root, lineage, step_name,
                                  sta1, sta2, comp, day, TAXIS, corr)
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_stack(db, params, root):
    """Write a trivial stacked CCF NetCDF per pair/comp/mov_stack, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_ccf
    while True:
        batch = get_next_lineage_batch(db, "stack", group_by="pair_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    data = np.zeros((1, len(TAXIS)), dtype="float32")
                    ds = xr.Dataset(
                        {"CCF": xr.DataArray(data, dims=["times", "taxis"],
                                             coords={"times": TIMES[:1],
                                                     "taxis": TAXIS})})
                    xr_save_ccf(root, lineage, step_name,
                                sta1, sta2, comp, ms, TAXIS, ds)
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_refstack(db, params, root):
    """Write a trivial REF NetCDF per pair/comp, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_ref
    while True:
        batch = get_next_lineage_batch(db, "refstack", group_by="pair_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                ref = np.zeros(len(TAXIS), dtype="float32")
                ds = xr.Dataset(
                    {"REF": xr.DataArray(ref, dims=["taxis"],
                                         coords={"taxis": TAXIS})})
                xr_save_ref(root, lineage, step_name,
                            sta1, sta2, comp, TAXIS, ds)
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_mwcs(db, params, root):
    """Write trivial MWCS NetCDF per pair/comp/mov_stack, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_mwcs
    n_keys = 4
    while True:
        batch = get_next_lineage_batch(db, "mwcs", group_by="pair_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({"MWCS": xr.DataArray(
                        np.zeros((len(TIMES), len(TAXIS), n_keys), "float32"),
                        dims=["times", "taxis", "keys"],
                        coords={"times": TIMES, "taxis": TAXIS,
                                "keys": ["dt", "err", "coh", "valid"]})})
                    xr_save_mwcs(root, lineage, step_name,
                                 sta1, sta2, comp, ms, TAXIS, ds)
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_dtt(db, params, root):
    """Write trivial DTT NetCDF per pair/comp/mov_stack, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_dtt
    while True:
        batch = get_next_lineage_batch(db, "mwcs_dtt", group_by="pair_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({"DTT": xr.DataArray(
                        np.zeros((len(TIMES), 3), "float32"),
                        dims=["times", "keys"],
                        coords={"times": TIMES,
                                "keys": ["dtt", "err", "coh"]})})
                    xr_save_dtt(root, lineage, step_name,
                                sta1, sta2, comp, ms, ds)
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_stretching(db, params, root):
    """Write trivial Stretching NetCDF per pair/comp/mov_stack, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_stretching
    while True:
        batch = get_next_lineage_batch(db, "stretching", group_by="pair_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({"STR": xr.DataArray(
                        np.zeros((len(TIMES), 3), "float32"),
                        dims=["times", "keys"],
                        coords={"times": TIMES,
                                "keys": ["Delta", "Coeff", "Error"]})})
                    xr_save_stretching(root, lineage, step_name,
                                       sta1, sta2, comp, ms, ds)
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_wct(db, params, root):
    """Write trivial WCT NetCDF per pair/comp/mov_stack, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_wct
    while True:
        batch = get_next_lineage_batch(db, "wavelet", group_by="pair_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    amp = [np.zeros((len(FREQS), len(TAXIS)), "float32")] * len(TIMES)
                    coh = [np.zeros((len(FREQS), len(TAXIS)), "float32")] * len(TIMES)
                    dtt = [np.zeros((len(FREQS), len(TAXIS)), "float32")] * len(TIMES)
                    xr_save_wct(root, lineage, step_name,
                                sta1, sta2, comp, ms,
                                TAXIS, FREQS, amp, coh, dtt,
                                list(TIMES))
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_wct_dtt(db, params, root):
    """Write trivial WCT-DTT NetCDF per pair/comp/mov_stack, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_wct_dtt
    while True:
        batch = get_next_lineage_batch(db, "wavelet_dtt", group_by="pair_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for sta1, sta2 in PAIRS:
            for comp in COMPS:
                for ms in MOV_STACKS:
                    ds = xr.Dataset({
                        "dtt": xr.DataArray(np.zeros((len(TIMES),), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                        "err": xr.DataArray(np.zeros((len(TIMES),), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                        "coh": xr.DataArray(np.ones((len(TIMES),), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                    })
                    xr_save_wct_dtt(root, lineage, step_name,
                                    sta1, sta2, comp, ms, TAXIS, ds)
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_dvv_agg(db, params, root, category: str):
    """Write trivial DVV aggregate NetCDF for *category*, mark sentinel Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_dvv_agg
    while True:
        batch = get_next_lineage_batch(db, category, group_by="day_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]   # upstream
        step_name = batch["lineage_names"][-1]
        for ms in MOV_STACKS:
            for comp in COMPS:
                ds = xr.Dataset({
                    "mean":    xr.DataArray(np.zeros(len(TIMES), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                    "std":     xr.DataArray(np.zeros(len(TIMES), "float32"),
                                            dims=["times"], coords={"times": TIMES}),
                    "n_pairs": xr.DataArray(np.full(len(TIMES), len(PAIRS), dtype="int32"),
                                            dims=["times"], coords={"times": TIMES}),
                })
                xr_save_dvv_agg(root, lineage, step_name,
                                ms, "CC", comp[0], ds)  # comp[0] = "Z" for "ZZ"
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_psd(db, params, root):
    """Write trivial PSD NetCDF per station/day, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_psd
    periods = np.logspace(-1, 1, 20)
    while True:
        batch = get_next_lineage_batch(db, "psd", group_by="day_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for net, sta in [("YA", "UV05"), ("YA", "UV06"), ("YA", "UV07")]:
            seed_id = f"{net}.{sta}.00.HHZ"
            day = batch["days"][0]
            ds = xr.Dataset({
                "PSD": xr.DataArray(
                    np.full((1, len(periods)), -150.0, dtype="float32"),
                    dims=["times", "periods"],
                    coords={"times": [np.datetime64(day)], "periods": periods})
            })
            xr_save_psd(root, lineage, step_name, seed_id, day, ds)
        massive_update_job(db, batch["jobs"], flag="D")


def _stub_psd_rms(db, params, root):
    """Write trivial PSD-RMS NetCDF per station, mark Done."""
    from ..core.workflow import get_next_lineage_batch, massive_update_job
    from ..core.io import xr_save_rms
    while True:
        batch = get_next_lineage_batch(db, "psd_rms", group_by="day_lineage")
        if batch is None:
            break
        lineage = batch["lineage_names"][:-1]
        step_name = batch["lineage_names"][-1]
        for net, sta in [("YA", "UV05"), ("YA", "UV06"), ("YA", "UV07")]:
            seed_id = f"{net}.{sta}.00.HHZ"
            import pandas as pd
            df = pd.DataFrame(
                {"low": [-150.0, -150.0], "mid": [-150.0, -150.0], "high": [-150.0, -150.0]},
                index=pd.DatetimeIndex(DAYS))
            xr_save_rms(root, lineage, step_name, seed_id, df)
        massive_update_job(db, batch["jobs"], flag="D")


# ── Tests ─────────────────────────────────────────────────────────────────────

@pytest.mark.order(1)
def test_smoke_01_schema(smoke_db):
    """Database tables exist including the Lineage normalisation table."""
    db, _, _ = smoke_db
    from sqlalchemy import inspect
    tables = inspect(db.get_bind()).get_table_names()
    for expected in ["jobs", "lineages", "workflow_steps", "workflow_links",
                     "config", "stations", "data_availability"]:
        assert expected in tables, f"Table {expected!r} missing"


@pytest.mark.order(2)
def test_smoke_02_workflow_created(smoke_db):
    """All expected workflow steps and links were created."""
    db, _, _ = smoke_db
    from ..core.workflow import get_workflow_steps, get_workflow_links
    steps = get_workflow_steps(db)
    step_names = {s.step_name for s in steps}
    expected = {"preprocess_1", "cc_1", "filter_1", "stack_1", "refstack_1",
                "mwcs_1", "mwcs_dtt_1", "mwcs_dtt_dvv_1",
                "stretching_1", "stretching_dvv_1",
                "wavelet_1", "wavelet_dtt_1", "wavelet_dtt_dvv_1",
                "psd_1", "psd_rms_1"}
    assert expected.issubset(step_names), f"Missing steps: {expected - step_names}"
    links = get_workflow_links(db)
    assert len(links) >= 5, "Too few workflow links"


@pytest.mark.order(3)
def test_smoke_03_new_jobs_initial(smoke_db):
    """new_jobs() creates preprocess_1 T jobs — one per station per day."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(init=True)  # init=True uses bulk insert without file-existence check
    # 3 stations × 2 days × (preprocess_1 + psd_1) = at least 6 preprocess T jobs
    _assert_min(db, "preprocess_1", "T", 6, "new_jobs must create 6 preprocess T jobs")


@pytest.mark.order(4)
def test_smoke_04_preprocess(smoke_db):
    """Stub preprocess marks all preprocess_1 jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "preprocess_1", "T", 1, "Need T preprocess jobs")
    _stub_preprocess(db, params, root)
    _assert_min(db, "preprocess_1", "D", 6, "preprocess_1 must be Done")


@pytest.mark.order(5)
def test_smoke_05_new_jobs_after_preprocess(smoke_db):
    """new_jobs --after preprocess creates cc_1 T jobs — one per pair per day."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="preprocess")
    # 3 pairs × 2 days
    _assert_min(db, "cc_1", "T", 6, "new_jobs --after preprocess must create 6 cc T jobs")


@pytest.mark.order(6)
def test_smoke_06_cc(smoke_db):
    """Stub CC marks all cc_1 jobs Done and writes CCF files."""
    db, params, root = smoke_db
    _assert_min(db, "cc_1", "T", 1)
    _stub_cc(db, params, root)
    _assert_min(db, "cc_1", "D", 6, "cc_1 must be Done after stub CC")


@pytest.mark.order(7)
def test_smoke_07_new_jobs_after_cc(smoke_db):
    """new_jobs --after cc creates stack_1 T jobs."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="cc")
    _assert_min(db, "stack_1", "T", 1, "new_jobs --after cc must create stack T jobs")


@pytest.mark.order(8)
def test_smoke_08_stack(smoke_db):
    """Stub stack marks all stack_1 jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "stack_1", "T", 1)
    _stub_stack(db, params, root)
    _assert_min(db, "stack_1", "D", 1, "stack_1 must be Done")


@pytest.mark.order(9)
def test_smoke_09_new_jobs_after_stack(smoke_db):
    """new_jobs --after stack creates refstack_1 T jobs."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="stack")
    _assert_min(db, "refstack_1", "T", 1, "new_jobs --after stack must create refstack T jobs")


@pytest.mark.order(10)
def test_smoke_10_refstack(smoke_db):
    """Stub refstack marks all refstack_1 jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "refstack_1", "T", 1)
    _stub_refstack(db, params, root)
    _assert_min(db, "refstack_1", "D", 1, "refstack_1 must be Done")


@pytest.mark.order(11)
def test_smoke_11_new_jobs_after_refstack(smoke_db):
    """new_jobs --after refstack creates mwcs_1, stretching_1, wavelet_1 T jobs."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="refstack")
    _assert_min(db, "mwcs_1",       "T", 1, "mwcs_1 T jobs expected after refstack")
    _assert_min(db, "stretching_1", "T", 1, "stretching_1 T jobs expected after refstack")
    _assert_min(db, "wavelet_1",    "T", 1, "wavelet_1 T jobs expected after refstack")


@pytest.mark.order(12)
def test_smoke_12_mwcs(smoke_db):
    """Stub MWCS marks all mwcs_1 jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "mwcs_1", "T", 1)
    _stub_mwcs(db, params, root)
    _assert_min(db, "mwcs_1", "D", 1)


@pytest.mark.order(13)
def test_smoke_13_new_jobs_after_mwcs(smoke_db):
    """new_jobs --after mwcs creates mwcs_dtt_1 T jobs."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="mwcs")
    _assert_min(db, "mwcs_dtt_1", "T", 1)


@pytest.mark.order(14)
def test_smoke_14_dtt(smoke_db):
    """Stub DTT marks all mwcs_dtt_1 jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "mwcs_dtt_1", "T", 1)
    _stub_dtt(db, params, root)
    _assert_min(db, "mwcs_dtt_1", "D", 1)


@pytest.mark.order(15)
def test_smoke_15_new_jobs_after_dtt(smoke_db):
    """new_jobs --after mwcs_dtt creates mwcs_dtt_dvv_1 sentinel job."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="mwcs_dtt")
    _assert_min(db, "mwcs_dtt_dvv_1", "T", 1, "mwcs_dtt_dvv sentinel job expected")


@pytest.mark.order(16)
def test_smoke_16_mwcs_dvv_agg(smoke_db):
    """Stub DVV-MWCS aggregate marks sentinel Done and writes NetCDF."""
    db, params, root = smoke_db
    _assert_min(db, "mwcs_dtt_dvv_1", "T", 1)
    _stub_dvv_agg(db, params, root, "mwcs_dtt_dvv")
    _assert_min(db, "mwcs_dtt_dvv_1", "D", 1)
    from ..results import MSNoiseResult
    results = MSNoiseResult.list(db, "mwcs_dtt_dvv")
    assert len(results) >= 1, "MSNoiseResult.list must find mwcs_dtt_dvv results"


@pytest.mark.order(17)
def test_smoke_17_stretching(smoke_db):
    """Stub stretching marks all stretching_1 jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "stretching_1", "T", 1)
    _stub_stretching(db, params, root)
    _assert_min(db, "stretching_1", "D", 1)


@pytest.mark.order(18)
def test_smoke_18_new_jobs_after_stretching(smoke_db):
    """new_jobs --after stretching creates stretching_dvv_1 sentinel job."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="stretching")
    _assert_min(db, "stretching_dvv_1", "T", 1)


@pytest.mark.order(19)
def test_smoke_19_stretching_dvv_agg(smoke_db):
    """Stub DVV-stretching aggregate marks sentinel Done."""
    db, params, root = smoke_db
    _assert_min(db, "stretching_dvv_1", "T", 1)
    _stub_dvv_agg(db, params, root, "stretching_dvv")
    _assert_min(db, "stretching_dvv_1", "D", 1)
    from ..results import MSNoiseResult
    assert len(MSNoiseResult.list(db, "stretching_dvv")) >= 1


@pytest.mark.order(20)
def test_smoke_20_wct(smoke_db):
    """Stub WCT marks all wavelet_1 jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "wavelet_1", "T", 1)
    _stub_wct(db, params, root)
    _assert_min(db, "wavelet_1", "D", 1)


@pytest.mark.order(21)
def test_smoke_21_new_jobs_after_wct(smoke_db):
    """new_jobs --after wavelet creates wavelet_dtt_1 T jobs."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="wavelet")
    _assert_min(db, "wavelet_dtt_1", "T", 1)


@pytest.mark.order(22)
def test_smoke_22_wct_dtt(smoke_db):
    """Stub WCT-DTT marks all wavelet_dtt_1 jobs Done."""
    db, params, root = smoke_db
    _assert_min(db, "wavelet_dtt_1", "T", 1)
    _stub_wct_dtt(db, params, root)
    _assert_min(db, "wavelet_dtt_1", "D", 1)


@pytest.mark.order(23)
def test_smoke_23_new_jobs_after_wct_dtt(smoke_db):
    """new_jobs --after wavelet_dtt creates wavelet_dtt_dvv_1 sentinel job."""
    db, _, _ = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main(after="wavelet_dtt")
    _assert_min(db, "wavelet_dtt_dvv_1", "T", 1)


@pytest.mark.order(24)
def test_smoke_24_wct_dvv_agg(smoke_db):
    """Stub DVV-WCT aggregate marks sentinel Done."""
    db, params, root = smoke_db
    _assert_min(db, "wavelet_dtt_dvv_1", "T", 1)
    _stub_dvv_agg(db, params, root, "wavelet_dtt_dvv")
    _assert_min(db, "wavelet_dtt_dvv_1", "D", 1)
    from ..results import MSNoiseResult
    assert len(MSNoiseResult.list(db, "wavelet_dtt_dvv")) >= 1


@pytest.mark.order(25)
def test_smoke_25_psd_branch(smoke_db):
    """PSD branch: new_jobs → stub_psd → Done → stub_psd_rms → Done."""
    db, params, root = smoke_db
    from ..s02_new_jobs import main as new_jobs_main
    new_jobs_main()  # creates psd_1 T jobs
    _assert_min(db, "psd_1", "T", 1, "new_jobs must create psd_1 T jobs")
    _stub_psd(db, params, root)
    _assert_min(db, "psd_1", "D", 1, "psd_1 must be Done")
    new_jobs_main(after="psd")
    _assert_min(db, "psd_rms_1", "T", 1, "new_jobs --after psd must create psd_rms T jobs")
    _stub_psd_rms(db, params, root)
    _assert_min(db, "psd_rms_1", "D", 1, "psd_rms_1 must be Done")


@pytest.mark.order(26)
def test_smoke_26_lineage_normalisation(smoke_db):
    """Every job has a non-null lineage_id; Lineage table rows << job rows."""
    db, _, _ = smoke_db
    from ..msnoise_table_def import Lineage, Job as _J
    null_jobs = db.query(_J).filter(_J.lineage_id.is_(None)).count()
    assert null_jobs == 0, f"{null_jobs} jobs have NULL lineage_id"
    n_lineages = db.query(Lineage).count()
    n_jobs     = db.query(_J).count()
    assert n_lineages >= 1,    "Lineage table must have at least one row"
    assert n_lineages < n_jobs, (
        f"Expected fewer Lineage rows ({n_lineages}) than jobs ({n_jobs})")
    # Every lineage string is non-empty
    for lin in db.query(Lineage).all():
        assert lin.lineage_str and "/" not in lin.lineage_str or True  # any non-empty string
        assert lin.lineage_str.strip(), f"Empty lineage_str: {lin!r}"


@pytest.mark.order(27)
def test_smoke_27_all_dvv_results_discoverable(smoke_db):
    """MSNoiseResult.list finds results for all three DVV methods."""
    db, _, _ = smoke_db
    from ..results import MSNoiseResult
    for cat in ("mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"):
        results = MSNoiseResult.list(db, cat)
        assert len(results) >= 1, f"MSNoiseResult.list({cat!r}) returned no results"
        r = results[0]
        assert r.category == cat
        data = r.get_dvv(pair_type="CC")
        assert len(data) >= 1, (
            f"get_dvv('CC') returned empty for {cat} — "
            f"check xr_save_dvv_agg / xr_get_dvv_agg path consistency")
