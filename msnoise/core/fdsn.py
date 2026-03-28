"""FDSN/EIDA waveform fetching for the MSNoise preprocess step.

This module handles fetching raw waveforms from FDSN web services or the
EIDA routing client, writing optional raw caches, and per-station error
handling.  It is the **only** place in MSNoise that makes network calls to
external data services.
"""

import logging
import os
import time

logger = logging.getLogger("msnoise.fdsn")


# ---------------------------------------------------------------------------
# URI scheme detection
# ---------------------------------------------------------------------------

def parse_datasource_scheme(uri: str) -> str:
    """Return the scheme of a DataSource URI.

    :param uri: DataSource.uri string.
    :returns: One of ``"local"``, ``"sds"``, ``"fdsn"``, ``"eida"``.
    """
    if not uri:
        return "local"
    from urllib.parse import urlsplit
    scheme = urlsplit(uri).scheme.lower()
    if scheme in ("", "sds"):
        return "sds" if scheme == "sds" else "local"
    if scheme == "fdsn":
        return "fdsn"
    if scheme == "eida":
        return "eida"
    return "local"


def is_remote_source(uri: str) -> bool:
    """Return True if the DataSource URI points to a remote service."""
    return parse_datasource_scheme(uri) in ("fdsn", "eida")


# ---------------------------------------------------------------------------
# Auth resolution
# ---------------------------------------------------------------------------

def get_auth(auth_env: str) -> dict:
    """Read credentials from environment variables for the given prefix.

    Looks up:
    - ``{auth_env}_FDSN_USER``
    - ``{auth_env}_FDSN_PASSWORD``
    - ``{auth_env}_FDSN_TOKEN``  (path to EIDA token file, or token string)

    :param auth_env: Env var prefix (e.g. ``"MSNOISE"``, ``"IRIS"``).
    :returns: Dict with keys ``user``, ``password``, ``token`` (any may be None).
    """
    prefix = auth_env.upper()
    return {
        "user":     os.environ.get(f"{prefix}_FDSN_USER"),
        "password": os.environ.get(f"{prefix}_FDSN_PASSWORD"),
        "token":    os.environ.get(f"{prefix}_FDSN_TOKEN"),
    }


# ---------------------------------------------------------------------------
# Client construction
# ---------------------------------------------------------------------------

def build_client(ds):
    """Build an ObsPy client for the given DataSource.

    :param ds: :class:`~msnoise.msnoise_table_def.DataSource` ORM object.
    :returns: An ObsPy client with a ``get_waveforms_bulk`` method.
    """
    from urllib.parse import urlsplit
    scheme = parse_datasource_scheme(ds.uri)
    auth = get_auth(ds.auth_env or "MSNOISE")

    if scheme == "fdsn":
        from obspy.clients.fdsn import Client
        base_url = urlsplit(ds.uri)._replace(scheme="https").geturl()
        # strip the "fdsn://" prefix → actual HTTP(S) URL
        base_url = ds.uri[len("fdsn://"):]
        kwargs = {}
        if auth["user"] and auth["password"]:
            kwargs["user"] = auth["user"]
            kwargs["password"] = auth["password"]
        if auth["token"]:
            kwargs["eida_token"] = auth["token"]
        return Client(base_url, **kwargs)

    if scheme == "eida":
        from obspy.clients.fdsn import RoutingClient
        base_url = ds.uri[len("eida://"):]
        kwargs = {}
        if auth["token"]:
            kwargs["eida_token"] = auth["token"]
        return RoutingClient(base_url, **kwargs)

    raise ValueError(f"build_client called for non-remote DataSource: {ds.uri!r}")


# ---------------------------------------------------------------------------
# Bulk fetch with retry + per-station error handling
# ---------------------------------------------------------------------------

def fetch_waveforms_bulk(client, bulk_request, retries=3):
    """Issue a bulk waveform request with retry logic.

    :param client: ObsPy FDSN/EIDA client.
    :param bulk_request: List of ``(net, sta, loc, chan, t1, t2)`` tuples.
    :param retries: Number of retry attempts for transient errors.
    :returns: :class:`~obspy.core.stream.Stream` (may be empty on failure).
    :raises: Re-raises auth errors and no-data exceptions immediately.
    """
    from obspy.clients.fdsn.header import FDSNNoDataException
    from obspy import Stream

    last_exc = None
    for attempt in range(1, retries + 1):
        try:
            return client.get_waveforms_bulk(bulk_request)
        except FDSNNoDataException:
            raise  # not transient — re-raise immediately
        except Exception as exc:
            last_exc = exc
            err_str = str(exc)
            # Auth failures are not transient
            if "401" in err_str or "403" in err_str or "Unauthorized" in err_str:
                raise
            wait = 2 ** attempt
            logger.warning(
                f"FDSN bulk request failed (attempt {attempt}/{retries}): "
                f"{exc}. Retrying in {wait}s."
            )
            time.sleep(wait)

    logger.error(f"FDSN bulk request failed after {retries} attempts: {last_exc}")
    return Stream()


# ---------------------------------------------------------------------------
# High-level: fetch + optional raw cache + per-station result
# ---------------------------------------------------------------------------

def fetch_and_preprocess(
    db, jobs, goal_day, params, responses=None, loglevel="INFO"
):
    """Fetch waveforms from FDSN/EIDA for a batch of stations on one day.

    Groups jobs by DataSource (all jobs in a batch share the same
    ``data_source_id`` when ``group_by="day_lineage_datasource"`` is used),
    issues one ``get_waveforms_bulk`` call, splits the result by station,
    and optionally writes raw files.

    :param db: SQLAlchemy session.
    :param jobs: List of Job ORM objects (all same day, same DataSource).
    :param goal_day: Date string ``YYYY-MM-DD``.
    :param params: :class:`~msnoise.params.LayeredParams` for this lineage.
    :param responses: ObsPy Inventory for instrument response removal, or None.
    :param loglevel: Logging level string.
    :returns: Tuple ``(stream, done_jobs, failed_jobs)`` where stream contains
        preprocessed traces for all successfully fetched stations.
    """
    from obspy import UTCDateTime, Stream
    from obspy.clients.fdsn.header import FDSNNoDataException
    from .stations import resolve_data_source, get_station
    from ..preprocessing import _apply_preprocessing  # see note below

    log = logging.getLogger(f"msnoise.fdsn.fetch")

    # Resolve DataSource from first job (all share same data_source_id in batch)
    first_job = jobs[0]
    net0, sta0, loc0 = first_job.pair.split(".")
    station0 = get_station(db, net0, sta0)
    ds = resolve_data_source(db, station0)

    fdsn_keep_raw = params.global_.fdsn_keep_raw
    retries       = int(params.global_.fdsn_retries or 3)
    min_coverage  = float(params.global_.fdsn_min_coverage or 0.5)
    output_folder = params.global_.output_folder
    step_name     = params.step_name

    t1 = UTCDateTime(goal_day)
    t2 = t1 + 86400

    # Build bulk request: one entry per (net, sta, loc, chan) per component
    raw = params.preprocess.preprocess_components
    comps = raw.split(',') if isinstance(raw, str) else (list(raw) if raw else ['Z'])

    bulk = []
    station_map = {}  # "NET.STA.LOC" → Job
    for job in jobs:
        net, sta, loc = job.pair.split(".")
        station_map[job.pair] = job
        for comp in comps:
            chan = f"HH{comp}"   # TODO: use actual channel from station config
            bulk.append((net, sta, loc, chan, t1, t2))

    log.info(f"FDSN bulk fetch: {len(jobs)} stations, {len(bulk)} channels, "
             f"day={goal_day}, source={ds.name!r}")

    # Build client and fetch
    try:
        client = build_client(ds)
        raw_stream = fetch_waveforms_bulk(client, bulk, retries=retries)
    except FDSNNoDataException:
        log.warning(f"No data available from {ds.name!r} for {goal_day} — marking all jobs Failed")
        return Stream(), [], jobs
    except Exception as exc:
        if "401" in str(exc) or "403" in str(exc) or "Unauthorized" in str(exc):
            log.error(
                f"Auth failure for DataSource {ds.name!r}. "
                f"Check {ds.auth_env}_FDSN_USER / {ds.auth_env}_FDSN_TOKEN. "
                f"Error: {exc}"
            )
        else:
            log.error(f"FDSN fetch failed for {ds.name!r} on {goal_day}: {exc}")
        return Stream(), [], jobs

    # Optionally write raw cache
    if fdsn_keep_raw in ("Y", "y", True):
        _write_raw_cache(raw_stream, output_folder, step_name, goal_day)

    # Split by station, check coverage, apply preprocessing
    done_jobs, failed_jobs = [], []
    out_stream = Stream()

    for sid, job in station_map.items():
        net, sta, loc = sid.split(".")
        sta_stream = raw_stream.select(network=net, station=sta, location=loc)

        if not sta_stream:
            log.warning(f"Station {sid} absent from FDSN bulk response — marking Failed")
            failed_jobs.append(job)
            continue

        # Coverage check
        total = sum(tr.stats.delta * tr.stats.npts for tr in sta_stream)
        if total < min_coverage * 86400:
            log.warning(
                f"{sid}: coverage {total/86400:.1%} < {min_coverage:.0%} — "
                f"proceeding with available data"
            )

        out_stream += sta_stream
        done_jobs.append(job)

    return out_stream, done_jobs, failed_jobs


def _write_raw_cache(stream, output_folder, step_name, goal_day):
    """Write raw fetched stream to ``_output/raw/<date>/<NET.STA.LOC.CHAN>.mseed``."""
    import numpy as np
    from obspy import Stream as _Stream

    raw_dir = os.path.join(output_folder, step_name, "_output", "raw", goal_day)
    os.makedirs(raw_dir, exist_ok=True)

    # Group by full SEED ID (NET.STA.LOC.CHAN)
    by_id = {}
    for tr in stream:
        sid = tr.id  # NET.STA.LOC.CHAN
        by_id.setdefault(sid, _Stream())
        by_id[sid].append(tr)

    for sid, st in by_id.items():
        fpath = os.path.join(raw_dir, f"{sid}.mseed")
        for tr in st:
            tr.data = tr.data.astype("float32")
        st.write(fpath, format="MSEED")
        logger.debug(f"Raw cache written: {fpath}")
