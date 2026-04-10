"""MSNoise xarray I/O for all result types (CCF, MWCS, DTT, STR, WCT, DVV, PSD)."""
import glob
import logging
import os
import sys

import numpy as np
import pandas as pd
import xarray as xr

# ── NetCDF encoding presets ───────────────────────────────────────────────────
#
# All float data variables are stored as float32 with zlib level-4 compression.
# float32 provides ~7 significant decimal digits — far exceeding measurement
# noise in any ambient-noise application — while reducing file size ~4–5×
# compared to uncompressed float64.
#
# Times (datetime64) and string/integer coordinates use their natural dtype.

_ZLIB_LEVEL = 4
_F32_OPTS   = {"dtype": "float32", "zlib": True, "complevel": _ZLIB_LEVEL}


def _f32_encoding(ds):
    """Build a float32+zlib encoding dict for all float data variables in *ds*."""
    enc = {}
    for var in ds.data_vars:
        da = ds[var]
        if da.dtype.kind == "f":           # float32 or float64
            enc[var] = _F32_OPTS.copy()
        elif da.dtype.kind in ("i", "u"):  # integers (e.g. n_pairs)
            enc[var] = {"zlib": True, "complevel": _ZLIB_LEVEL}
    return enc


def _trim(data, dttname, limits=0.1):
    """
    Trimmed mean and standard deviation calculation.

    :param data: DataFrame containing the data.
    :param dttname: Name of the column used for grouping.
    :param limits: Trimming limits (default is 0.1).
    :return: Tuple containing the trimmed mean and trimmed standard deviation.
    """
    from scipy.stats.mstats import trimmed_mean, trimmed_std
    grouped = data[dttname].groupby(level=0)
    if limits == 0:
        g = grouped.mean()
        h = grouped.std()
    else:
        g = grouped.apply(trimmed_mean, limits=limits)
        h = grouped.apply(trimmed_std, limits=limits)
    return g, h

# ============================================================
# ── Core helpers ───────────────────────────────────────────


def _xr_create_or_open(fn, taxis=[], name="CCF", lazy=False):
    if os.path.isfile(fn):
        # Write paths (save functions) must use load_dataset so the file handle
        # is fully released before _xr_save_and_close writes back to the same
        # file.  Read-only paths pass lazy=True to get a memory-mapped dataset
        # that avoids loading all data into RAM upfront.
        if lazy:
            return xr.open_dataset(fn)
        ds = xr.load_dataset(fn)
        return ds
    times = pd.date_range("2000-01-01", freq="h", periods=0)
    if name == "CCF":
        data = np.random.random((len(times), len(taxis)))
        dr = xr.DataArray(data, coords=[times, taxis], dims=["times", "taxis"])
    elif name == "REF":
        data = np.random.random(len(taxis))
        dr = xr.DataArray(data, coords=[taxis], dims=["taxis"])
    elif name == "MWCS":
        keys = ["M", "EM", "MCOH"]
        data = np.random.random((len(times), len(taxis), len(keys)))
        dr = xr.DataArray(data, coords=[times, taxis, keys],
                          dims=["times", "taxis", "keys"])
    elif name == "STR":
        keys = ["Delta", "Coeff", "Error"]
        data = np.random.random((len(times), len(keys)))
        dr = xr.DataArray(data, coords=[times, keys],
                          dims=["times", "keys"])
    elif name == "DTT":
        keys = ["m", "em", "a", "ea", "m0", "em0"]
        data = np.random.random((len(times), len(keys)))
        dr = xr.DataArray(data, coords=[times, keys],
                          dims=["times", "keys"])
    elif name == "DVV":
        level0 = ["m", "em", "a", "ea", "m0", "em0"]
        level1 = ['10%', '25%', '5%', '50%', '75%', '90%', '95%', 'count', 'max', 'mean',
                  'min', 'std', 'trimmed_mean', 'trimmed_std', 'weighted_mean', 'weighted_std']
        data = np.random.random((len(times), len(level0), len(level1)))
        dr = xr.DataArray(data, coords=[times, level0, level1],
                          dims=["times", "level0", "level1"])
    elif name == "WCT":
        dvv_data = np.random.random((len(times), len(taxis)))
        err_data = np.random.random((len(times), len(taxis)))
        coh_data = np.random.random((len(times), len(taxis)))

        dvv_da = xr.DataArray(dvv_data, coords=[times, taxis], dims=['times', 'frequency'])
        err_da = xr.DataArray(err_data, coords=[times, taxis], dims=['times', 'frequency'])
        coh_da = xr.DataArray(coh_data, coords=[times, taxis], dims=['times', 'frequency'])

        # Combine into a Dataset
        ds = xr.Dataset({'DTT': dvv_da, 'ERR': err_da, 'COH': coh_da})
        return ds

    else:
        logging.error("Not implemented, name=%s invalid." % name)
        sys.exit(1)
    dr.name = name
    return dr.to_dataset()



def _xr_insert_or_update(dataset, new):
    tt = new.merge(dataset, compat='override', combine_attrs="drop_conflicts", join='outer')
    return tt.combine_first(dataset)



def _xr_save_and_close(dataset, fn, encoding=None):
    """Write *dataset* to *fn* with float32+zlib encoding by default."""
    if not os.path.isdir(os.path.split(fn)[0]):
        os.makedirs(os.path.split(fn)[0], exist_ok=True)
    enc = encoding if encoding is not None else _f32_encoding(dataset)
    dataset.to_netcdf(fn, mode="w", encoding=enc, engine="netcdf4")
    dataset.close()
    del dataset






def xr_save_ccf(root, lineage, step_name, station1, station2, components, mov_stack, taxis, new, overwrite=False):
    path = os.path.join(root, *lineage, step_name, "_output",
                        "%s_%s" % (mov_stack[0], mov_stack[1]), "%s" % components)
    fn = "%s_%s.nc" % (station1, station2)
    fullpath = os.path.join(path, fn)
    if overwrite:
        _xr_save_and_close(new, fullpath)
    else:
        dr = _xr_create_or_open(fullpath, taxis, name="CCF")
        dr = _xr_insert_or_update(dr, new)
        dr = dr.sortby("times")
        _xr_save_and_close(dr, fullpath)
        return dr



def xr_get_ccf(root, lineage, station1, station2, components, mov_stack, taxis, format="dataset"):
    """Load CCF results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataset"`` (default) returns a :class:`xarray.DataArray` (CCF variable).
        ``"dataframe"`` returns a :class:`~pandas.DataFrame` (legacy).
    """
    path = os.path.join(root, *lineage, "_output",
                        "%s_%s" % (mov_stack[0], mov_stack[1]), "%s" % components)
    fn = "%s_%s.nc" % (station1, station2)

    fullpath = os.path.join(path, fn)
    if not os.path.isfile(fullpath):
        raise FileNotFoundError(fullpath)
    data = _xr_create_or_open(fullpath, taxis, name="CCF", lazy=True)
    if format == "dataset":
        return data.CCF
    # ── DataFrame (legacy) ──────────────────────────────────────────────
    return data.CCF.to_dataframe().unstack().droplevel(0, axis=1)




def xr_save_ccf_daily(root, lineage, step_name, station1, station2, components, date, taxis, corr):
    """Save a single daily-stacked CCF as a per-day NetCDF file.

    One file per calendar day — safe for concurrent workers since each worker
    owns exactly one day.

    Path layout::

        <root>/<lineage>/<step_name>/_output/daily/<components>/<sta1>_<sta2>/<YYYY-MM-DD>.nc

    :param corr: 1-D numpy array of the stacked CCF.
    :param taxis: 1-D lag-time axis array.
    :param date: ``datetime.date`` or ISO string ``"YYYY-MM-DD"``.
    """
    path = os.path.join(root, *lineage, step_name, "_output",
                        "daily", components,
                        f"{station1}_{station2}")
    os.makedirs(path, exist_ok=True)
    fn = os.path.join(path, f"{date}.nc")
    da = xr.DataArray(corr, coords=[taxis], dims=["taxis"], name="CCF")
    _xr_save_and_close(da.to_dataset(name=da.name or "CCF"), fn)



def xr_get_ccf_daily(root, lineage, step_name, station1, station2, components, date):
    """Load a single daily-stacked CCF written by :func:`xr_save_ccf_daily`.

    :returns: :class:`xarray.DataArray` with dim ``taxis``.
    :raises FileNotFoundError: if the NetCDF file does not exist.
    """
    fn = os.path.join(root, *lineage, step_name, "_output",
                      "daily", components,
                      f"{station1}_{station2}", f"{date}.nc")
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    ds = xr.load_dataset(fn)
    return ds["CCF"]



def xr_save_ref(root, lineage, step_name, station1, station2, components, taxis, new, overwrite=False):
    path = os.path.join(root, *lineage, step_name, "_output",
                        "REF", "%s" % components)
    fn = "%s_%s.nc" % (station1, station2)
    fullpath = os.path.join(path, fn)
    if overwrite:
        _xr_save_and_close(new, fullpath)
    else:
        dr = _xr_create_or_open(fullpath, taxis, name="REF")
        dr = _xr_insert_or_update(dr, new)
        _xr_save_and_close(dr, fullpath)
        return dr



def xr_get_ref(root, lineage, station1, station2, components, taxis, ignore_network=False):
    path = os.path.join(root, *lineage, "_output",
                        "REF", "%s" % components)
    # If ignore_network is True, strip the network code from the station names
    if ignore_network:
        s1_parts = station1.split('.')
        s2_parts = station2.split('.')

        available_files = glob.glob(os.path.join(path, "*.%s.%s_*.%s.%s.nc" % (s1_parts[1],s1_parts[2], s2_parts[1], s1_parts[2])))

        if available_files:
            # Use the first available reference file
            fullpath = available_files[0]
        else:
            raise FileNotFoundError(f"No reference file found for station {s1_parts[1]} and {s2_parts[1]}")
    else:
        fn = "%s_%s.nc" % (station1, station2)

        fullpath = os.path.join(path, fn)
        if not os.path.isfile(fullpath):
            # logging.error("FILE DOES NOT EXIST: %s, skipping" % fullpath)
            raise FileNotFoundError(fullpath)
    data = _xr_create_or_open(fullpath, taxis, name="REF", lazy=True)
    return data


def xr_load_ccf_for_stack(root, lineage_names, station1, station2, components, dates):
    """Load per-day CCF data for stacking — the higher-level replacement for ``get_results_all``.

    Reads the daily CCF NetCDF files written by :func:`~msnoise.s03_compute_no_rotation`
    for the requested station pair, component and date list.  Tries ``keep_all``
    (per-window, ``_output/all/``) first; falls back to ``keep_days``
    (daily stacks, ``_output/daily/``) if the former are absent.

    Path layout (written by s03)::

        <root> / *cc_lineage / filter_step / _output / all|daily / <comp> / <sta1>_<sta2> / <date>.nc

    :param root: Output folder (``params.global_.output_folder``).
    :param lineage_names: Full lineage name list including the filter step,
        e.g. ``["preprocess_1", "cc_1", "filter_1", "stack_1"]``.
    :param station1: First station SEED id ``NET.STA.LOC``.
    :param station2: Second station SEED id ``NET.STA.LOC``.
    :param components: Component pair string e.g. ``"ZZ"``.
    :param dates: Iterable of ``datetime.date`` or ISO date strings.
    :returns: :class:`xarray.Dataset` with variable ``CCF`` and dims
        ``(times, taxis)``, sorted by ``times``.  Empty Dataset if no data found.
    """
    # Derive cc_lineage and filter_step from lineage_names
    cc_idx = next(
        (i for i, n in enumerate(lineage_names) if n.startswith("cc_")), None
    )
    if cc_idx is None or cc_idx + 1 >= len(lineage_names):
        return xr.Dataset()
    cc_lineage = lineage_names[:cc_idx + 1]   # e.g. ['preprocess_1', 'cc_1']
    filter_step = lineage_names[cc_idx + 1]   # e.g. 'filter_1'

    das = []
    for date in dates:
        date_str = date if isinstance(date, str) else date.strftime('%Y-%m-%d')
        try:
            da = xr_get_ccf_all(root, cc_lineage, filter_step,
                                station1, station2, components, date_str)
            das.append(da)
        except FileNotFoundError:
            pass

    if not das:
        # Fallback: daily stacks from keep_days (xr_save_ccf_daily)
        for date in dates:
            date_str = date if isinstance(date, str) else date.strftime('%Y-%m-%d')
            try:
                da = xr_get_ccf_daily(root, cc_lineage, filter_step,
                                      station1, station2, components, date_str)
                t = pd.Timestamp(date_str)
                da = da.expand_dims({"times": [t]})
                das.append(da)
            except FileNotFoundError:
                pass

    if not das:
        return xr.Dataset()
    combined = xr.concat(das, dim="times").sortby("times")
    return combined.to_dataset(name="CCF")


# ── MWCS ────────────────────────────────────────────────────


def xr_save_mwcs(root, lineage, step_name, station1, station2, components, mov_stack, taxis, dataset):
    """Save MWCS results to a NetCDF file.

    :param dataset: :class:`xarray.Dataset` with a ``MWCS`` variable and
        dims ``(times, taxis, keys)``, as built by
        :mod:`~msnoise.s05_compute_mwcs`.
    """
    fn = os.path.join(root, *lineage, step_name, "_output",
                        "%s_%s" % (mov_stack[0], mov_stack[1]),
                       "%s" % components,
                       "%s_%s.nc" % (station1, station2))
    os.makedirs(os.path.split(fn)[0], exist_ok=True)
    dr = _xr_create_or_open(fn, taxis=dataset.coords.get("taxis", taxis), name="MWCS")
    rr = _xr_insert_or_update(dr, dataset)
    _xr_save_and_close(rr, fn)



def xr_get_mwcs(root, lineage, station1, station2, components, mov_stack, format="dataset"):
    """Load MWCS results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`
        with MultiIndex columns ``(keys, taxis)``.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset`.
    """
    fn = os.path.join(root, *lineage, "_output",
                        "%s_%s" % (mov_stack[0], mov_stack[1]),
                       "%s" % components,
                       "%s_%s.nc" % (station1, station2))
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    data = _xr_create_or_open(fn, name="MWCS", lazy=True)

    if format == "dataset":
        return data

    # ── DataFrame (legacy) ──────────────────────────────────────
    da = data.MWCS  # DataArray with dims (times, taxis, keys)
    # Build DataFrame directly — pandas-version-safe
    times_vals = da.coords["times"].values
    taxis_vals = da.coords["taxis"].values
    keys_vals  = da.coords["keys"].values
    n_t, n_tx, n_k = da.values.shape
    # Transpose to (times, keys, taxis) so MultiIndex is (keys, taxis)
    midx = pd.MultiIndex.from_product([keys_vals, taxis_vals], names=["keys", "taxis"])
    arr = da.values.transpose(0, 2, 1).reshape(n_t, n_k * n_tx)
    return pd.DataFrame(arr, index=pd.DatetimeIndex(times_vals), columns=midx)


# ── DTT ─────────────────────────────────────────────────────


def xr_save_dtt(root, lineage, step_name, station1, station2, components, mov_stack, dataset):
    """Save DTT results to a NetCDF file.

    :param dataset: :class:`xarray.Dataset` with a ``DTT`` variable and
        dims ``(times, keys)``, as built by :mod:`~msnoise.s06_compute_mwcs_dtt`.
    """
    fn = os.path.join(root, *lineage, step_name, "_output",
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components,
                      "%s_%s.nc" % (station1, station2))
    os.makedirs(os.path.split(fn)[0], exist_ok=True)
    dr = _xr_create_or_open(fn, taxis=[], name="DTT")
    rr = _xr_insert_or_update(dr, dataset)
    _xr_save_and_close(rr, fn)



def xr_get_dtt(root, lineage, station1, station2, components, mov_stack, format="dataset"):
    """Load DTT results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset`.
    """
    fn = os.path.join(root, *lineage, "_output",
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      "%s" % components,
                      "%s_%s.nc" % (station1, station2))

    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    dr = _xr_create_or_open(fn, taxis=[], name="DTT", lazy=True)

    if format == "dataset":
        return dr

    # ── DataFrame (legacy) ──────────────────────────────────────────────
    da = dr.DTT  # DataArray with dims (times, keys)
    return pd.DataFrame(
        da.values,
        index=pd.DatetimeIndex(da.coords["times"].values),
        columns=list(da.coords["keys"].values),
    )


# ── Stretching ──────────────────────────────────────────────


def xr_save_stretching(root, lineage, step_name, station1, station2,
                        components, mov_stack, dataset):
    """Save per-pair stretching results to a NetCDF file.

    Path layout::

        <root>/<lineage>/<step_name>/_output/<mov_stack[0]>_<mov_stack[1]>/<components>/<sta1>_<sta2>.nc

    :param dataset: :class:`xarray.Dataset` with a ``STR`` variable and
        dims ``(times, keys)`` where keys = ``['Delta', 'Coeff', 'Error']``,
        as built by :mod:`~msnoise.s10_stretching`.
    """
    fn = os.path.join(
        root, *lineage, step_name, "_output",
        "%s_%s" % (mov_stack[0], mov_stack[1]),
        components,
        "%s_%s.nc" % (station1, station2),
    )
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    dr = _xr_create_or_open(fn, taxis=[], name="STR")
    rr = _xr_insert_or_update(dr, dataset)
    _xr_save_and_close(rr, fn)



def _xr_get_stretching(root, lineage, station1, station2, components, mov_stack, format="dataset"):
    """Load per-pair stretching results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`
        with columns ``Delta``, ``Coeff``, ``Error``.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset`.

    :raises FileNotFoundError: if the NetCDF file does not exist.
    """
    fn = os.path.join(
        root, *lineage, "_output",
        "%s_%s" % (mov_stack[0], mov_stack[1]),
        components,
        "%s_%s.nc" % (station1, station2),
    )
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    dr = _xr_create_or_open(fn, taxis=[], name="STR", lazy=True)

    if format == "dataset":
        return dr

    # ── DataFrame (legacy) ──────────────────────────────────────────────
    da = dr.STR  # DataArray with dims (times, keys)
    return pd.DataFrame(
        da.values,
        index=pd.DatetimeIndex(da.coords["times"].values),
        columns=list(da.coords["keys"].values),
    )


# ── WCT ─────────────────────────────────────────────────────


def xr_save_wct(root, lineage, step_name, station1, station2, components, mov_stack, taxis, freqs, WXamp_list, WXcoh_list, WXdt_list, dates_list):
    """
    Save WCT results into an xarray Dataset and store it as a NetCDF file.

    Parameters:
    - root: str, Root output folder path
    - lineage: list, Lineage path components
    - step_name: str, Step name for output path
    - station1, station2: str, Station pair
    - components: str, Seismic component (e.g., ZZ)
    - mov_stack: tuple, Moving stack window (e.g., ('1d', '1d'))
    - taxis, freqs: np.array, Time axis and frequency axis
    - WXamp_list, WXcoh_list, WXdt_list: list of np.array, WCT outputs
    - dates_list: list of datetime, Timestamps for each WCT calculation
    """

    # Convert lists to xarray DataArrays (all WCT outputs are real-valued;
    # cast explicitly in case intermediate complex dtype was inherited)
    WXamp_da = xr.DataArray(
        data=np.array(WXamp_list).real,
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="WXamp"
    )

    Wcoh_da = xr.DataArray(
        data=np.array(WXcoh_list).real,
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="Wcoh"
    )

    WXdt_da = xr.DataArray(
        data=np.array(WXdt_list).real,
        dims=["times", "freqs", "taxis"],
        coords={"times": dates_list, "freqs": freqs, "taxis": taxis},
        name="WXdt"
    )

    # Combine into an xarray Dataset
    ds = xr.Dataset({"WXamp": WXamp_da, "Wcoh": Wcoh_da, "WXdt": WXdt_da})

    # Define output directory
    fn = os.path.join(root, *lineage, step_name, "_output",
                      "%s_%s" % (mov_stack[0], mov_stack[1]),
                      components, f"{station1}_{station2}.nc")

    os.makedirs(os.path.dirname(fn), exist_ok=True)

    # Save to NetCDF
    _xr_save_and_close(ds, fn)


    # Cleanup memory
    del ds, WXamp_da, Wcoh_da, WXdt_da



def xr_load_wct(root, lineage, station1, station2, components, mov_stack):
    """
    Load WCT results from an xarray Dataset stored in a NetCDF file.

    Parameters:
    - root: str, Root output folder path
    - lineage: list, Lineage path components
    - station1, station2: str, Station pair
    - components: str, Seismic component (e.g., ZZ)
    - mov_stack: tuple, Moving stack window (e.g., ('1d', '1d'))

    Returns:
    - ds: xarray.Dataset containing the WCT data (WXamp, Wcoh, WXdt)
    """

    # Construct the file path
    fn = os.path.join(root, *lineage, "_output",
                      f"{mov_stack[0]}_{mov_stack[1]}", components,
                      f"{station1}_{station2}.nc")

    # Check if the file exists
    if not os.path.exists(fn):
        raise FileNotFoundError(f"File not found: {fn}")

    # Load and return the dataset lazily — caller reads but never writes back
    ds = xr.open_dataset(fn)
    return ds



def xr_save_wct_dtt(root, lineage, step_name, station1, station2, components, mov_stack, taxis, dataset):
    """Save WCT-DTT results to a NetCDF file.

    :param dataset: :class:`xarray.Dataset` with variables ``DTT``, ``ERR``,
        ``COH`` and dims ``(times, frequency)``, as built by
        :mod:`~msnoise.s09_compute_wct_dtt`.
    """
    fn = os.path.join(root, *lineage, step_name, "_output",
                      f"{mov_stack[0]}_{mov_stack[1]}", components,
                      f"{station1}_{station2}.nc")
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    existing_ds = _xr_create_or_open(fn, name="WCT")
    updated_ds = _xr_insert_or_update(existing_ds, dataset)
    _xr_save_and_close(updated_ds, fn)
    logging.debug(f"Saved WCT DTT data to {fn}")



def xr_get_wct_dtt(root, lineage, station1, station2, components, mov_stack):
    """Load per-pair WCT dt/t results from a NetCDF file.

    Returns an :class:`xarray.Dataset` with variables ``DTT``, ``ERR``,
    ``COH`` each of shape ``(times, frequency)``, as written by
    :func:`xr_save_wct_dtt`.

    :raises FileNotFoundError: if the NetCDF file does not exist.
    :rtype: :class:`xarray.Dataset`
    """
    fn = os.path.join(
        root, *lineage, "_output",
        "%s_%s" % (mov_stack[0], mov_stack[1]),
        components,
        "%s_%s.nc" % (station1, station2),
    )
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    return xr.open_dataset(fn)



def xr_save_ccf_all(root, lineage, step_name, station1, station2,
                    components, date, window_times, taxis, corrs):
    """Save all windowed CCFs for one day as a single NetCDF file.

    Replaces the legacy HDF5-based ``export_allcorr``.  One file per calendar
    day — safe for concurrent workers since each worker owns exactly one day.

    Path layout::

        <root>/<lineage>/<step_name>/_output/all/<components>/<sta1>_<sta2>/<date>.nc

    :param window_times: 1-D array of window start times (datetime strings or
        datetime objects) — the ``times`` dimension.
    :param taxis: 1-D lag-time axis array — the ``taxis`` dimension.
    :param corrs: 2-D numpy array ``(n_windows, n_taxis)`` of CCF data.
    """
    path = os.path.join(root, *lineage, step_name, "_output",
                        "all", components,
                        f"{station1}_{station2}")
    os.makedirs(path, exist_ok=True)
    fn = os.path.join(path, f"{date}.nc")
    times = pd.to_datetime(window_times)
    da = xr.DataArray(
        corrs,
        coords=[times, taxis],
        dims=["times", "taxis"],
        name="CCF",
    )
    _xr_save_and_close(da.to_dataset(name=da.name or "CCF"), fn)



def xr_get_ccf_all(root, lineage, step_name, station1, station2,
                   components, date):
    """Load all windowed CCFs for one day written by :func:`xr_save_ccf_all`.

    :returns: :class:`xarray.DataArray` with dims ``(times, taxis)``.
    :raises FileNotFoundError: if the NetCDF file does not exist.
    """
    fn = os.path.join(root, *lineage, step_name, "_output",
                      "all", components,
                      f"{station1}_{station2}", f"{date}.nc")
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    ds = xr.load_dataset(fn)
    return ds["CCF"]




def save_daily_ccf(root, lineage, step_name, station1, station2,
                   components, date, corr, taxis):
    """Save a daily-stacked CCF to NetCDF using :func:`xr_save_ccf_daily`.

    Replaces the legacy :func:`add_corr` / ``_export_mseed`` / ``_export_sac``
    pipeline.  Output lives at::

        <root>/<lineage>/<step_name>/_output/daily/<components>/<sta1>_<sta2>/<date>.nc

    :param root: Output folder (``params.global_.output_folder``).
    :param lineage: List of step-name strings for the lineage path.
    :param step_name: Current step name (e.g. ``"cc_1"``).
    :param station1: First station SEED id ``NET.STA.LOC``.
    :param station2: Second station SEED id ``NET.STA.LOC``.
    :param components: Component pair string e.g. ``"ZZ"``.
    :param date: Calendar date (``datetime.date`` or ISO string).
    :param corr: 1-D numpy array, the stacked CCF.
    :param taxis: 1-D lag-time axis array.
    """
    xr_save_ccf_daily(root, lineage, step_name, station1, station2,
                      components, date, taxis, corr)





def _classify_pair_type(sta1: str, sta2: str, component: str = "") -> str:
    """Classify a station pair as CC, SC, or AC.

    Classification is based on the SEED ``NET.STA.LOC`` ids and the component
    pair string (e.g. ``"ZZ"``, ``"ZN"``):

    - **AC** (autocorrelation): ``sta1 == sta2`` *and* both component letters
      are identical (e.g. ``"ZZ"``, ``"EE"``).  Same instrument correlating
      with itself.
    - **SC** (single-channel cross-component): ``sta1 == sta2`` *and*
      component letters differ (e.g. ``"ZN"``).  Different channels at the
      same location.
    - **CC** (cross-correlation): ``sta1 != sta2``.  Different physical
      locations regardless of component.

    :param sta1: First station SEED id ``NET.STA.LOC``.
    :param sta2: Second station SEED id ``NET.STA.LOC``.
    :param component: Component-pair string (e.g. ``"ZZ"``, ``"ZN"``).
        If empty or length < 2, falls back to comparing only the SEED ids.
    :returns: ``"AC"``, ``"SC"``, or ``"CC"``.
    """
    if sta1 != sta2:
        return "CC"
    # Same NET.STA.LOC — distinguish by component
    if len(component) >= 2 and component[0] == component[-1]:
        return "AC"
    return "SC"



def _dvv_column_spec(parent_category: str, pair_type: str, params) -> tuple:
    """Return ``(dv_col, err_col, quality_col)`` for the given parent category.

    For ``mwcs_dtt``, the columns depend on the pair type and user config
    (smart defaults: CC→free-intercept m/em, SC/AC→forced-zero m0/em0).
    For ``stretching`` and ``wct_dtt`` the columns are fixed.

    :param parent_category: ``"mwcs_dtt"``, ``"stretching"``, or ``"wavelet_dtt"``.
    :param pair_type: ``"CC"``, ``"SC"``, ``"AC"``, or ``"ALL"``.
    :param params: Params object from :func:`get_params`.
    :returns: Tuple ``(dv_col, err_col, quality_col)`` where any of the
        quality columns may be ``None`` if not available.
    """
    if parent_category == "mwcs_dtt":
        pt = pair_type if pair_type in ("CC", "SC", "AC") else "CC"
        pt_lower = pt.lower()
        dv_col  = getattr(params.category_layer, f"dvv_{pt_lower}_value",  "m" if pt == "CC" else "m0")
        err_col = getattr(params.category_layer, f"dvv_{pt_lower}_error", "em" if pt == "CC" else "em0")
        return dv_col, err_col, "mcoh"
    elif parent_category == "stretching":
        return "Delta", "Error", "Coeff"
    elif parent_category == "wavelet_dtt":
        # WCT stores DTT/ERR/COH per frequency; scalar extracted separately
        return "DTT", "ERR", "COH"
    else:
        raise ValueError(f"Unknown parent_category: {parent_category!r}")



def _freq_average_wct(ds, freqmin: float, freqmax: float,
                       quality_min: float, freq_agg: str = "mean"):
    """Collapse the ``(times, frequency)`` WCT-DTT Dataset to scalar time series.

    Selects the frequency band ``[freqmin, freqmax]``, masks cells where
    ``COH < quality_min``, then collapses the frequency axis using *freq_agg*.

    :param ds: :class:`xarray.Dataset` with variables ``DTT``, ``ERR``, ``COH``
        and dims ``(times, frequency)``.
    :param freqmin: Lower frequency bound (Hz).
    :param freqmax: Upper frequency bound (Hz).
    :param quality_min: Minimum coherence; cells below are set to NaN.
    :param freq_agg: ``"mean"`` or ``"median"`` over the frequency axis.
    :returns: Tuple ``(dv_da, err_da)`` — two 1-D :class:`xarray.DataArray`
        with dim ``times``.
    :raises ValueError: if no frequencies fall within [freqmin, freqmax].
    """
    freqs = ds.coords["frequency"].values
    mask_freq = (freqs >= freqmin) & (freqs <= freqmax)
    if not mask_freq.any():
        raise ValueError(
            f"No frequencies in [{freqmin}, {freqmax}] Hz in WCT-DTT output. "
            f"Available: {freqs.min():.3f}–{freqs.max():.3f} Hz"
        )
    dtt_sel = ds["DTT"].isel(frequency=mask_freq)
    err_sel = ds["ERR"].isel(frequency=mask_freq)
    coh_sel = ds["COH"].isel(frequency=mask_freq)

    # Mask low-coherence cells
    bad = coh_sel < quality_min
    dtt_sel = dtt_sel.where(~bad)
    err_sel = err_sel.where(~bad)

    # Collapse frequency axis
    if freq_agg == "median":
        dv_da  = dtt_sel.median(dim="frequency")
        err_da = err_sel.median(dim="frequency")
    else:
        dv_da  = dtt_sel.mean(dim="frequency")
        err_da = err_sel.mean(dim="frequency")

    return dv_da, err_da



def xr_save_dvv_agg(root, lineage, step_name, mov_stack,
                    pair_type: str, component: str, dataset):
    """Save a DVV aggregate result to a NetCDF file.

    Path layout::

        <root>/<lineage>/<step_name>/_output/<mov_stack>/dvv_<pair_type>_<component>.nc

    :param dataset: :class:`xarray.Dataset` with dim ``times`` and stat
        variables as built by :func:`aggregate_dvv_pairs`.
    """
    ms_str = "%s_%s" % (mov_stack[0], mov_stack[1])
    fn = os.path.join(root, *lineage, step_name, "_output",
                      ms_str, f"dvv_{pair_type}_{component}.nc")
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    _xr_save_and_close(dataset, fn)



def xr_get_dvv_agg(root, lineage, step_name, mov_stack,
                   pair_type: str, component: str, format: str = "dataset"):
    """Load a DVV aggregate result from a NetCDF file.

    :param format: ``"dataset"`` (default) or ``"dataframe"``.
    :raises FileNotFoundError: if the file does not exist.
    """
    ms_str = "%s_%s" % (mov_stack[0], mov_stack[1])
    fn = os.path.join(root, *lineage, step_name, "_output",
                      ms_str, f"dvv_{pair_type}_{component}.nc")
    if not os.path.isfile(fn):
        raise FileNotFoundError(fn)
    ds = xr.open_dataset(fn)

    if format == "dataframe":
        # materialise before returning as DataFrame — lazy → eager here is fine
        return ds.load().to_dataframe()
    return ds



def aggregate_dvv_pairs(root, parent_lineage, parent_step_name,
                        parent_category: str, mov_stack, component: str,
                        pair_type: str, pairs, params) -> xr.Dataset:
    """Aggregate per-pair DTT/STR/WCT-DTT results into network-level dv/v statistics.

    Reads all per-pair output files for the given ``(mov_stack, component)``
    combination, extracts the appropriate dv/v and error columns per pair type,
    applies quality filtering, then computes across-pair statistics at each
    time step.

    :param root: Output folder root.
    :param parent_lineage: Lineage name list up to and including the parent
        DTT step (e.g. ``["preprocess_1", ..., "mwcs_dtt_1"]``).
    :param parent_step_name: Step name of the parent DTT step.
    :param parent_category: ``"mwcs_dtt"``, ``"stretching"``, or ``"wavelet_dtt"``.
    :param mov_stack: Tuple ``(window, step)`` e.g. ``("1D", "1D")``.
    :param component: Component string e.g. ``"ZZ"``.
    :param pair_type: ``"CC"``, ``"SC"``, ``"AC"``, or ``"ALL"``.
    :param pairs: Iterable of ``(sta1, sta2)`` SEED-id string tuples.
    :param params: Params object from :func:`get_params`.
    :returns: :class:`xarray.Dataset` with dim ``times`` and stat variables.
    :raises ValueError: if no data files are found for the given combination.
    """
    dv_col, err_col, quality_col = _dvv_column_spec(
        parent_category, pair_type, params)

    # All xr_get_* functions expect *lineage to already end with the step name,
    # e.g. [..., "mwcs_dtt_1"].  parent_lineage = lineage_names_upstream of the
    # DVV job, which strips the DVV step and ends with the DTT step name —
    # exactly what the getters need.  No adjustment required.
    dtt_lineage = list(parent_lineage)

    quality_min    = params.category_layer.dvv_quality_min
    do_weighted    = str(params.category_layer.dvv_weighted_mean).upper() == "Y"
    do_trimmed     = str(params.category_layer.dvv_trimmed_mean).upper() == "Y"
    trim_sigma     = params.category_layer.dvv_trim_limit
    out_percent    = str(params.category_layer.dvv_output_percent).upper() == "Y"
    percentile_str = params.category_layer.dvv_percentiles
    percentiles    = [float(p) for p in str(percentile_str).split(",") if p.strip()]

    # WCT-specific params
    if parent_category == "wavelet_dtt":
        wct_freqmin  = params.category_layer.dvv_freqmin
        wct_freqmax  = params.category_layer.dvv_freqmax
        wct_freq_agg = str(params.category_layer.dvv_freq_agg)

    # ── 1. Collect per-pair 1-D time series ─────────────────────────────
    pair_dvv  = []   # list of (times_array, dv_array, err_array)

    for sta1, sta2 in pairs:
        # Filter by pair_type
        pt = _classify_pair_type(sta1, sta2, component)
        if pair_type != "ALL" and pt != pair_type:
            continue

        try:
            if parent_category == "mwcs_dtt":
                ds = xr_get_dtt(root, dtt_lineage, sta1, sta2,
                                 component, mov_stack, format="dataset")
                # dv/v = -1 * dt/t estimated by compute_mwcs_dtt
                da_dv  = -1 * ds["DTT"].sel(keys=dv_col)
                da_err = ds["DTT"].sel(keys=err_col)
                if quality_col in ds["DTT"].coords["keys"].values:
                    da_q = ds["DTT"].sel(keys=quality_col)
                    bad = (da_q < quality_min).values  # numpy bool; scalar 'keys' coord on da_dv/da_q
                    da_dv  = da_dv.copy(data=np.where(bad, np.nan, da_dv.values))
                    da_err = da_err.copy(data=np.where(bad, np.nan, da_err.values))

            elif parent_category == "stretching":
                ds = _xr_get_stretching(root, dtt_lineage, sta1, sta2,
                                         component, mov_stack, format="dataset")
                # dv/v is "stretching factor - 1" because in msnoise the REF is stretched, not the current CCF.
                da_dv  = ds["STR"].sel(keys=dv_col) - 1
                da_err = ds["STR"].sel(keys=err_col)
                if quality_col:
                    da_q = ds["STR"].sel(keys=quality_col)
                    bad = (da_q < quality_min).values  # numpy bool; scalar 'keys' coord on da_dv/da_q
                    da_dv  = da_dv.copy(data=np.where(bad, np.nan, da_dv.values))
                    da_err = da_err.copy(data=np.where(bad, np.nan, da_err.values))

            elif parent_category == "wavelet_dtt":
                ds = xr_get_wct_dtt(root, dtt_lineage, sta1, sta2,
                                     component, mov_stack)
                da_dv, da_err = _freq_average_wct(
                    ds, wct_freqmin, wct_freqmax, quality_min, wct_freq_agg)
                # wavelet dtt is returning dt/t -> *=-1 for dv/v
                da_dv *= -1

        except FileNotFoundError as e:
            logging.getLogger("msnoise").warning(
                f"DTT file not found for {sta1}:{sta2}/{component} "
                f"(lineage={dtt_lineage}, mov_stack={mov_stack}): {e}")
            continue
        except Exception as e:
            import traceback
            logging.getLogger("msnoise").warning(
                f"Error reading {sta1}:{sta2}/{component}: {e}\n"
                + traceback.format_exc()
            )
            continue

        if out_percent:
            da_dv = da_dv * 100.0

        times = da_dv.coords["times"].values
        dv    = da_dv.values.astype(float)
        err   = da_err.values.astype(float)
        pair_dvv.append((times, dv, err))

    if not pair_dvv:
        raise ValueError(
            f"No data found for parent={parent_category} "
            f"step={parent_step_name} mov_stack={mov_stack} "
            f"component={component} pair_type={pair_type} "
            f"in folder={dtt_lineage}"
        )

    # ── 2. Build (pairs × times) arrays on a common time axis ───────────
    all_times_sorted = np.array(sorted({t for times, _, _ in pair_dvv
                                         for t in times}),
                                dtype="datetime64[ns]")
    n_t = len(all_times_sorted)
    n_p = len(pair_dvv)

    dv_mat  = np.full((n_p, n_t), np.nan)
    err_mat = np.full((n_p, n_t), np.nan)

    time_idx = {t: i for i, t in enumerate(all_times_sorted)}
    for p, (times, dv, err) in enumerate(pair_dvv):
        for ti, t in enumerate(times):
            key = np.datetime64(t, "ns")
            if key in time_idx:
                j = time_idx[key]
                dv_mat[p, j]  = dv[ti]
                err_mat[p, j] = err[ti]

    # ── 3. Compute statistics across pairs at each time step ─────────────
    # n_pairs: count of non-NaN pairs per time step
    n_pairs = np.sum(~np.isnan(dv_mat), axis=0).astype(float)

    mean_dv   = np.nanmean(dv_mat, axis=0)
    std_dv    = np.nanstd(dv_mat, axis=0, ddof=1)
    median_dv = np.nanmedian(dv_mat, axis=0)

    data_vars = {
        "mean":     ("times", mean_dv),
        "std":      ("times", std_dv),
        "median":   ("times", median_dv),
        "n_pairs":  ("times", n_pairs),
    }

    # Percentiles
    for pct in percentiles:
        key = f"q{int(pct):02d}"
        pct_vals = np.nanpercentile(dv_mat, pct, axis=0)
        data_vars[key] = ("times", pct_vals)

    # Weighted mean/std
    if do_weighted:
        w_mat = np.where(err_mat > 0, 1.0 / err_mat**2, 0.0)
        w_sum = np.nansum(w_mat, axis=0)
        w_sum_safe = np.where(w_sum > 0, w_sum, np.nan)
        wmean = np.nansum(w_mat * np.nan_to_num(dv_mat), axis=0) / w_sum_safe
        # weighted std
        wvar  = np.nansum(w_mat * (np.nan_to_num(dv_mat) - wmean[None, :])**2,
                           axis=0) / w_sum_safe
        wstd  = np.sqrt(wvar)
        data_vars["weighted_mean"] = ("times", wmean)
        data_vars["weighted_std"]  = ("times", wstd)

    # Trimmed mean/std (sigma-based: remove pairs > trim_sigma*std from mean)
    if do_trimmed:
        tmean = np.full(n_t, np.nan)
        tstd  = np.full(n_t, np.nan)
        for j in range(n_t):
            col = dv_mat[:, j]
            valid = col[~np.isnan(col)]
            if len(valid) < 2:
                tmean[j] = np.nan if len(valid) == 0 else valid[0]
                tstd[j]  = np.nan
                continue
            mu = np.mean(valid)
            sigma = np.std(valid, ddof=1)
            kept = valid[np.abs(valid - mu) <= trim_sigma * sigma]
            tmean[j] = np.mean(kept) if len(kept) > 0 else np.nan
            tstd[j]  = np.std(kept, ddof=1) if len(kept) > 1 else np.nan
        data_vars["trimmed_mean"] = ("times", tmean)
        data_vars["trimmed_std"]  = ("times", tstd)

    # ── 4. Build output Dataset ──────────────────────────────────────────
    ds_out = xr.Dataset(
        {k: xr.DataArray(v[1], dims=["times"],
                          coords={"times": all_times_sorted})
         for k, v in data_vars.items()},
        attrs={
            "parent_category": parent_category,
            "parent_step":     parent_step_name,
            "pair_type":       pair_type,
            "component":       component,
            "mov_stack":       "%s_%s" % (mov_stack[0], mov_stack[1]),
            "dvv_unit":        "percent" if out_percent else "fraction",
            "created":         str(np.datetime64("now")),
        }
    )
    return ds_out

# ============================================================


def psd_read_results(net, sta, loc, chan, datelist, format='PPSD', use_cache=True):
    from obspy.signal import PPSD
    if loc == "--":
        loc = ""
    fn = "%s.%s.%s.%s-%s_%s.npz" % (net, sta, loc, chan, datelist[0], datelist[-1])
    import tempfile
    fn = os.path.join(tempfile.gettempdir(), "MSNOISE-PSD", fn)
    if use_cache and os.path.isfile(fn):
        logging.debug("I found this cool file: %s" % fn)
        ppsd = PPSD.load_npz(fn)
    else:
        first = True
        ppsd = None
        for day in datelist:
            jday = int(day.strftime("%j"))
            toglob = os.path.join('PSD', 'NPZ', "%s" % day.year, net, sta,
                                  chan + ".D", "%s.%s.%s.%s.D.%s.%03i.npz" % (
                                  net, sta, loc, chan, day.year, jday))
            files = glob.glob(toglob)
            if not len(files):
                logging.error("No files found for %s.%s.%s.%s: %s" % (
                net, sta, loc, chan, day))
                continue
            file = files[0]
            if os.path.isfile(file):
                if first:
                    ppsd = PPSD.load_npz(file)
                    first = False
                else:
                    try:
                        ppsd.add_npz(file)
                    except Exception:
                        pass
    if not ppsd:
        return None
    if use_cache:
        if not os.path.isdir(os.path.split(fn)[0]):
            os.makedirs(os.path.split(fn)[0], exist_ok=True)
        ppsd.save_npz(fn[:-4])
    return ppsd



def psd_ppsd_to_dataframe(ppsd):
    from obspy import UTCDateTime
    ind_times = np.array(
        [UTCDateTime(t).datetime for t in ppsd.current_times_used])
    data = np.asarray(ppsd._binned_psds)
    return pd.DataFrame(data, index=ind_times, columns=ppsd.period_bin_centers)



def psd_ppsd_to_dataset(ppsd):
    """Convert an ObsPy PPSD object to an :class:`xarray.Dataset`.

    Builds the same data as :func:`psd_ppsd_to_dataframe` but returns an
    ``xr.Dataset`` with a ``PSD`` variable of dims ``(times, periods)`` —
    ready to pass directly to :func:`xr_save_psd` without a DataFrame
    round-trip.

    :param ppsd: :class:`~obspy.signal.spectral_estimation.PPSD` object.
    :returns: :class:`xarray.Dataset`.
    """
    from obspy import UTCDateTime
    times = np.array([UTCDateTime(t).datetime for t in ppsd.current_times_used])
    periods = np.asarray(ppsd.period_bin_centers, dtype=float)
    data = np.asarray(ppsd._binned_psds, dtype=float)
    da = xr.DataArray(
        data,
        coords=[times, periods],
        dims=["times", "periods"],
        name="PSD",
    )
    return da.to_dataset()



def xr_save_psd(root, lineage, step_name, seed_id, day, dataset):
    """Save a daily PSD result to a NetCDF file.

    Path layout::

        <root>/<lineage>/<step_name>/_output/daily/<seed_id>/<YYYY-MM-DD>.nc

    :param dataset: :class:`xarray.Dataset` with a ``PSD`` variable and
        dims ``(times, periods)``, as built by :func:`psd_ppsd_to_dataset`.
    """
    day_str = day if isinstance(day, str) else day.strftime("%Y-%m-%d")
    fn = os.path.join(
        root, *lineage, step_name, "_output", "daily",
        seed_id, f"{day_str}.nc",
    )
    os.makedirs(os.path.dirname(fn), exist_ok=True)
    _xr_save_and_close(dataset, fn)



def xr_load_psd(root, lineage, step_name, seed_id, day, format="dataset"):
    """Load a daily PSD NetCDF written by :func:`xr_save_psd`.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`
        or ``None`` if file not found.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset` or ``None``.
    """
    day_str = day if isinstance(day, str) else day.strftime("%Y-%m-%d")
    fn = os.path.join(
        root, *lineage, step_name, "_output", "daily",
        seed_id, f"{day_str}.nc",
    )
    if not os.path.isfile(fn):
        return None
    ds = xr.open_dataset(fn)

    if format == "dataset":
        return ds

    # ── DataFrame (legacy) ──────────────────────────────────────
    da = ds.PSD
    df = pd.DataFrame(
        da.values,
        index=pd.DatetimeIndex(da.coords["times"].values),
        columns=da.coords["periods"].values.astype(float),
    )
    ds.close()
    return df


def xr_save_rms(root, lineage, step_name, seed_id, dataframe):
    """Save per-station PSD RMS results to a NetCDF file.

    Accepts either a :class:`~pandas.DataFrame` (legacy, index = DatetimeIndex,
    columns = band labels) or an :class:`xarray.Dataset` with a ``RMS``
    variable and dims ``(times, bands)``.
    """
    fn = os.path.join(root, *lineage, step_name, "_output", seed_id, "RMS.nc")
    os.makedirs(os.path.dirname(fn), exist_ok=True)

    # ── Dataset path (new, xarray-native) ──────────────────────────────
    if isinstance(dataframe, xr.Dataset):
        if os.path.isfile(fn):
            existing_ds = xr.load_dataset(fn)
            merged = _xr_insert_or_update(existing_ds, dataframe)
            _xr_save_and_close(merged, fn)
        else:
            _xr_save_and_close(dataframe, fn)
        return

    # ── DataFrame path (legacy compat) ─────────────────────────────────
    bands = list(dataframe.columns.astype(str))
    times = pd.DatetimeIndex(dataframe.index)

    if os.path.isfile(fn):
        existing = xr_load_rms(root, lineage, step_name, seed_id,
                               format="dataframe")
        if existing is not None:
            existing = existing[~existing.index.isin(times)]
            dataframe = pd.concat([existing, dataframe]).sort_index()
            bands = list(dataframe.columns.astype(str))
            times = pd.DatetimeIndex(dataframe.index)

    da = xr.DataArray(
        dataframe.values.astype(float),
        coords=[times, bands],
        dims=["times", "bands"],
        name="RMS",
    )
    _xr_save_and_close(da.to_dataset(name=da.name or "data"), fn)



def xr_load_rms(root, lineage, step_name, seed_id, format="dataset"):
    """Load per-station PSD RMS results from a NetCDF file.

    Parameters
    ----------
    format : str
        ``"dataframe"`` returns a :class:`~pandas.DataFrame`
        or ``None`` if file not found.
        ``"dataset"`` (default) returns an :class:`xarray.Dataset` or ``None``.
    """
    fn = os.path.join(root, *lineage, step_name, "_output", seed_id, "RMS.nc")
    if not os.path.isfile(fn):
        return None
    ds = xr.open_dataset(fn)

    if format == "dataset":
        return ds

    # ── DataFrame (legacy) ──────────────────────────────────────
    df = pd.DataFrame(
        ds.RMS.values,
        index=pd.DatetimeIndex(ds.RMS.coords["times"].values),
        columns=list(ds.RMS.coords["bands"].values),
    )
    ds.close()
    return df
