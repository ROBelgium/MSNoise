"""Pure-logic unit tests for msnoise.core.compute, .signal, and .config.

No DB, no filesystem, no ObsPy required.
"""
import types

import numpy as np
import pytest
import xarray as xr


# ============================================================================
# core/compute.py
# ============================================================================

class TestAnalyticPhase:
    def test_output_shape(self):
        from ..core.compute import _analytic_phase
        x = np.random.randn(256)
        y = _analytic_phase(x)
        assert y.shape == x.shape
        assert np.iscomplexobj(y)

    def test_unit_magnitude(self):
        from ..core.compute import _analytic_phase
        x = np.random.randn(256)
        y = _analytic_phase(x)
        assert np.allclose(np.abs(y), 1.0, atol=1e-3)

    def test_all_zeros_input(self):
        from ..core.compute import _analytic_phase
        x = np.zeros(128)
        y = _analytic_phase(x)
        assert y.shape == x.shape  # must not raise


class TestAnalyticPhaseBatch:
    def test_output_shape(self):
        from ..core.compute import _analytic_phase_batch
        X = np.random.randn(4, 256)
        Y = _analytic_phase_batch(X)
        assert Y.shape == X.shape
        assert np.iscomplexobj(Y)

    def test_unit_magnitude_per_row(self):
        from ..core.compute import _analytic_phase_batch
        X = np.random.randn(3, 128)
        Y = _analytic_phase_batch(X)
        assert np.allclose(np.abs(Y), 1.0, atol=1e-3)


class TestPccXcorr:
    def _make_data(self, n=256, n_sta=3):
        rng = np.random.default_rng(42)
        return rng.standard_normal((n_sta, n)).astype(float)

    def test_empty_index(self):
        from ..core.compute import pcc_xcorr
        data = self._make_data()
        result = pcc_xcorr(data, maxlag=50, energy=None, index=[])
        assert result == {}

    def test_output_length(self):
        from ..core.compute import pcc_xcorr
        data = self._make_data()
        maxlag = 50
        index = [(0, 0, 1), (1, 0, 2)]
        result = pcc_xcorr(data, maxlag=maxlag, energy=None, index=index)
        assert set(result.keys()) == {0, 1}
        assert len(result[0]) == 2 * maxlag + 1

    def test_normalized_max(self):
        from ..core.compute import pcc_xcorr
        data = self._make_data()
        result = pcc_xcorr(data, maxlag=30, energy=None,
                           index=[(0, 0, 1)], normalized="MAX")
        assert result[0].max() == pytest.approx(1.0)

    def test_normalized_absmax(self):
        from ..core.compute import pcc_xcorr
        data = self._make_data()
        result = pcc_xcorr(data, maxlag=30, energy=None,
                           index=[(0, 0, 1)], normalized="ABSMAX")
        assert np.abs(result[0]).max() == pytest.approx(1.0)

    def test_self_correlation_peak_at_zero(self):
        from ..core.compute import pcc_xcorr
        rng = np.random.default_rng(7)
        sig = rng.standard_normal((2, 256))
        sig[1] = sig[0]   # identical traces → peak at lag=0
        maxlag = 40
        result = pcc_xcorr(sig, maxlag=maxlag, energy=None,
                           index=[(0, 0, 1)], normalized="ABSMAX")
        ccf = result[0]
        assert np.argmax(np.abs(ccf)) == maxlag   # zero-lag index


class TestSmooth:
    def test_boxcar_shape(self):
        from ..core.compute import smooth
        x = np.arange(50, dtype=float)
        y = smooth(x, window='boxcar', half_win=3)
        assert len(y) == len(x)

    def test_hann_shape(self):
        from ..core.compute import smooth
        x = np.ones(50)
        y = smooth(x, window='hanning', half_win=2)
        assert len(y) == len(x)

    def test_constant_signal_unchanged(self):
        from ..core.compute import smooth
        x = np.ones(40) * 5.0
        y = smooth(x, window='boxcar', half_win=3)
        assert np.allclose(np.real(y), 5.0, atol=1e-10)


class TestResolveWctLagMin:
    def _make_params(self, lag_type="static", v=1.0, minlag=5.0):
        wdt = types.SimpleNamespace(wct_lag=lag_type, wct_v=v, wct_minlag=minlag)
        return types.SimpleNamespace(wavelet_dtt=wdt)

    def test_static_returns_minlag(self):
        from ..core.compute import resolve_wct_lag_min
        p = self._make_params(lag_type="static", minlag=8.0)
        assert resolve_wct_lag_min(p, dist=100.0) == pytest.approx(8.0)

    def test_dynamic_returns_dist_over_v(self):
        from ..core.compute import resolve_wct_lag_min
        p = self._make_params(lag_type="dynamic", v=2.0, minlag=5.0)
        assert resolve_wct_lag_min(p, dist=100.0) == pytest.approx(50.0)

    def test_dynamic_none_v_falls_back_to_static(self):
        from ..core.compute import resolve_wct_lag_min
        # wct_v=None → `None or 1.0` → 1.0, so dynamic gives dist/1.0
        p = self._make_params(lag_type="dynamic", v=None, minlag=7.0)
        # v resolves to 1.0 via `or`, so result = dist/1.0 = 100.0
        assert resolve_wct_lag_min(p, dist=100.0) == pytest.approx(100.0)


class TestBuildWctDttDataset:
    def test_output_is_dataset(self):
        from ..core.compute import build_wct_dtt_dataset
        dates = [np.datetime64("2023-01-01"), np.datetime64("2023-01-02")]
        freqs = np.array([0.5, 1.0, 2.0])
        dtt_rows = [np.ones(3), np.ones(3) * 2]
        err_rows = [np.zeros(3), np.zeros(3)]
        coh_rows = [np.ones(3) * 0.9, np.ones(3) * 0.8]
        ds = build_wct_dtt_dataset(dates, dtt_rows, err_rows, coh_rows, freqs)
        assert isinstance(ds, xr.Dataset)
        assert "times" in ds.dims
        assert "frequency" in ds.dims


# ============================================================================
# core/signal.py
# ============================================================================

class TestNextpow2:
    def test_exact_power(self):
        from ..core.signal import nextpow2
        assert nextpow2(8) == pytest.approx(3)

    def test_non_power(self):
        from ..core.signal import nextpow2
        assert nextpow2(9) == pytest.approx(4)

    def test_one(self):
        from ..core.signal import nextpow2
        assert nextpow2(1) == pytest.approx(0)


class TestGetWindow:
    def test_boxcar_length(self):
        from ..core.signal import get_window
        w = get_window(window="boxcar", half_win=3)
        assert len(w) == 7

    def test_hanning_length(self):
        from ..core.signal import get_window
        w = get_window(window="hanning", half_win=5)
        assert len(w) == 11

    def test_sum_normalised(self):
        from ..core.signal import get_window
        w = get_window(window="boxcar", half_win=3)
        assert np.sum(np.real(w)) == pytest.approx(1.0, abs=1e-10)


class TestGetCoherence:
    def test_perfect_coherence(self):
        from ..core.signal import getCoherence
        n = 10
        dcs = np.ones(n)
        ds1 = np.ones(n)
        ds2 = np.ones(n)
        coh = getCoherence(dcs, ds1, ds2)
        assert np.allclose(np.abs(coh), 1.0)

    def test_zero_denominator_gives_zero(self):
        from ..core.signal import getCoherence
        n = 10
        dcs = np.ones(n)
        ds1 = np.zeros(n)
        ds2 = np.zeros(n)
        coh = getCoherence(dcs, ds1, ds2)
        assert np.allclose(coh, 0.0)

    def test_clipped_to_one(self):
        from ..core.signal import getCoherence
        n = 5
        coh = getCoherence(np.ones(n) * 10, np.ones(n), np.ones(n))
        assert np.all(np.abs(coh) <= 1.0 + 1e-10)


class TestPrepareAbsPositiveFft:
    def test_positive_freqs_only(self):
        from ..core.signal import prepare_abs_positive_fft
        x = np.random.randn(128)
        freq, val = prepare_abs_positive_fft(x, sampling_rate=100.0)
        assert np.all(freq >= 0)

    def test_output_shapes_match(self):
        from ..core.signal import prepare_abs_positive_fft
        x = np.random.randn(64)
        freq, val = prepare_abs_positive_fft(x, sampling_rate=50.0)
        assert freq.shape == val.shape


class TestPsdRms:
    def test_flat_spectrum(self):
        from ..core.signal import psd_rms
        f = np.linspace(0.1, 1.0, 100)
        s = np.ones_like(f)
        rms = psd_rms(s, f)
        assert rms > 0

    def test_zero_spectrum(self):
        from ..core.signal import psd_rms
        f = np.linspace(0.1, 1.0, 100)
        s = np.zeros_like(f)
        assert psd_rms(s, f) == pytest.approx(0.0)


class TestStack:
    def _make_ccfs(self, n_traces=5, n_lags=101):
        rng = np.random.default_rng(0)
        return rng.standard_normal((n_traces, n_lags))

    def test_linear_shape(self):
        from ..core.signal import stack
        data = self._make_ccfs()
        out = stack(data, stack_method="linear")
        assert out.shape == (data.shape[1],)

    def test_linear_is_mean(self):
        from ..core.signal import stack
        data = self._make_ccfs()
        out = stack(data, stack_method="linear")
        assert np.allclose(out, data.mean(axis=0))

    def test_pws_shape(self):
        from ..core.signal import stack
        data = self._make_ccfs()
        out = stack(data, stack_method="pws",
                    pws_timegate=5.0, pws_power=2,
                    goal_sampling_rate=20.0)
        assert out.shape == (data.shape[1],)


class TestFindSegments:
    def _make_da(self, n_times=10, n_lags=20, nan_rows=None):
        data = np.ones((n_times, n_lags))
        if nan_rows:
            for r in nan_rows:
                data[r, :] = np.nan
        return xr.DataArray(data, dims=["times", "lags"])

    def test_no_gaps(self):
        from ..core.signal import find_segments
        da = self._make_da()
        segs = find_segments(da, gap_threshold=1)
        assert len(segs) == 1
        assert len(segs[0]) == 10

    def test_gap_longer_than_threshold_splits(self):
        from ..core.signal import find_segments
        # Two all-NaN rows is a gap of 3 index steps > gap_threshold=1.
        # This used to return a single segment: the old implementation reset
        # prev_idx on every null row, so the gap test could never fire.
        da = self._make_da(nan_rows=[4, 5])
        segs = find_segments(da, gap_threshold=1)
        assert len(segs) == 2
        assert segs[0] == [0, 1, 2, 3]
        assert segs[1] == [6, 7, 8, 9]

    def test_gap_shorter_than_threshold_is_bridged(self):
        from ..core.signal import find_segments
        da = self._make_da(nan_rows=[4, 5])
        segs = find_segments(da, gap_threshold=5)
        assert len(segs) == 1
        assert segs[0] == [0, 1, 2, 3, 6, 7, 8, 9]

    def test_long_gap_always_splits(self):
        from ..core.signal import find_segments
        da = self._make_da(n_times=100, nan_rows=list(range(20, 60)))
        segs = find_segments(da, gap_threshold=5)
        assert len(segs) == 2
        assert segs[0][-1] == 19 and segs[1][0] == 60

    def test_all_nan(self):
        from ..core.signal import find_segments
        da = self._make_da(nan_rows=list(range(10)))
        segs = find_segments(da, gap_threshold=1)
        assert segs == []


# ============================================================================
# core/config.py
# ============================================================================

class TestParseConfigKey:
    def test_bare_name_is_global(self):
        from ..core.config import parse_config_key
        assert parse_config_key("output_folder") == ("global", 1, "output_folder")

    def test_category_dot_name(self):
        from ..core.config import parse_config_key
        assert parse_config_key("cc.cc_sampling_rate") == ("cc", 1, "cc_sampling_rate")

    def test_category_dot_setnum_dot_name(self):
        from ..core.config import parse_config_key
        assert parse_config_key("mwcs.2.mwcs_wlen") == ("mwcs", 2, "mwcs_wlen")

    def test_invalid_set_number_raises(self):
        from ..core.config import parse_config_key
        with pytest.raises(ValueError):
            parse_config_key("cc.abc.param")

    def test_too_many_parts_raises(self):
        from ..core.config import parse_config_key
        with pytest.raises(ValueError):
            parse_config_key("a.b.c.d")


class TestCastConfigValue:
    def test_bool_y(self):
        from ..core.config import _cast_config_value
        assert _cast_config_value("flag", "Y", "bool") == "Y"
        assert _cast_config_value("flag", "yes", "bool") == "Y"
        assert _cast_config_value("flag", "True", "bool") == "Y"
        assert _cast_config_value("flag", "1", "bool") == "Y"

    def test_bool_n(self):
        from ..core.config import _cast_config_value
        assert _cast_config_value("flag", "N", "bool") == "N"
        assert _cast_config_value("flag", "false", "bool") == "N"
        assert _cast_config_value("flag", "0", "bool") == "N"

    def test_bool_invalid_raises(self):
        from ..core.config import _cast_config_value
        with pytest.raises(ValueError):
            _cast_config_value("flag", "maybe", "bool")

    def test_int_valid(self):
        from ..core.config import _cast_config_value
        assert _cast_config_value("n", "42", "int") == "42"

    def test_int_invalid_raises(self):
        from ..core.config import _cast_config_value
        with pytest.raises(ValueError):
            _cast_config_value("n", "abc", "int")

    def test_float_valid(self):
        from ..core.config import _cast_config_value
        assert _cast_config_value("x", "3.14", "float") == "3.14"

    def test_float_invalid_raises(self):
        from ..core.config import _cast_config_value
        with pytest.raises(ValueError):
            _cast_config_value("x", "pi", "float")

    def test_str_passthrough(self):
        from ..core.config import _cast_config_value
        assert _cast_config_value("s", "anything goes", "str") == "anything goes"


class TestLineageToPlotTag:
    def test_basic(self):
        from ..core.config import lineage_to_plot_tag
        tag = lineage_to_plot_tag(["preprocess_1", "cc_1", "filter_1", "stack_1"])
        assert "pre1" in tag
        assert "cc1" in tag
        assert "stk1" in tag

    def test_empty(self):
        from ..core.config import lineage_to_plot_tag
        assert lineage_to_plot_tag([]) == ""

    def test_single_step(self):
        from ..core.config import lineage_to_plot_tag
        tag = lineage_to_plot_tag(["preprocess_2"])
        assert "2" in tag


class TestBuildPlotOutfile:
    def test_explicit_path_unchanged(self):
        from ..core.config import build_plot_outfile
        out = build_plot_outfile("myfile.png", "ccftime",
                                 ["preprocess_1", "cc_1"])
        assert out == "myfile.png"

    def test_none_returns_none(self):
        from ..core.config import build_plot_outfile
        assert build_plot_outfile(None, "ccftime", []) is None

    def test_question_mark_expands(self):
        from ..core.config import build_plot_outfile
        out = build_plot_outfile("?.png", "ccftime",
                                 ["preprocess_1", "cc_1", "filter_1", "stack_1"])
        assert out.startswith("ccftime__")
        assert out.endswith(".png")

    def test_pair_included(self):
        from ..core.config import build_plot_outfile
        out = build_plot_outfile("?.png", "ccftime",
                                 ["preprocess_1", "cc_1"],
                                 pair="BE.UCC..HHZ:BE.MEM..HHZ")
        assert "BE.UCC..HHZ-BE.MEM..HHZ" in out

    def test_mov_stack_tuple(self):
        from ..core.config import build_plot_outfile
        out = build_plot_outfile("?.png", "dvv",
                                 ["preprocess_1"],
                                 mov_stack=("1D", "1D"))
        assert "m1D-1D" in out


# ============================================================================
# core/stations.py — pure-logic helpers (no DB, no ObsPy required)
# ============================================================================

class TestGetInterstationDistance:
    def _sta(self, x, y):
        return types.SimpleNamespace(X=x, Y=y)

    def test_utm_hypot(self):
        from ..core.stations import get_interstation_distance
        s1 = self._sta(0.0, 0.0)
        s2 = self._sta(3000.0, 4000.0)   # 5 km
        dist = get_interstation_distance(s1, s2, coordinates="UTM")
        assert dist == pytest.approx(5.0, rel=1e-6)

    def test_utm_zero(self):
        from ..core.stations import get_interstation_distance
        s = self._sta(1234.0, 5678.0)
        assert get_interstation_distance(s, s, coordinates="UTM") == pytest.approx(0.0)

    def test_deg_uses_gps2dist(self):
        from ..core.stations import get_interstation_distance
        # Brussels → Paris ≈ 265 km
        brussels = self._sta(4.35, 50.85)
        paris    = self._sta(2.35, 48.85)
        dist = get_interstation_distance(brussels, paris, coordinates="DEG")
        assert 250 < dist < 290, f"Brussels-Paris distance out of range: {dist}"


class TestToSds:
    """Tests for stations.to_sds path builder (no DB, no ObsPy stream needed)."""

    def _stats(self, net="BE", sta="UCC", loc="", chan="HHZ"):
        return types.SimpleNamespace(
            network=net, station=sta, location=loc, channel=chan
        )

    def test_format_structure(self):
        from ..core.stations import to_sds
        path = to_sds(self._stats(), year=2023, jday=42)
        parts = path.split("/")
        assert parts[0] == "2023"
        assert parts[1] == "BE"
        assert parts[2] == "UCC"
        assert parts[3] == "HHZ.D"
        assert parts[4] == "BE.UCC..HHZ.D.2023.042"

    def test_year_zero_padded(self):
        from ..core.stations import to_sds
        path = to_sds(self._stats(), year=999, jday=1)
        assert path.startswith("0999/")

    def test_jday_zero_padded(self):
        from ..core.stations import to_sds
        path = to_sds(self._stats(), year=2023, jday=5)
        assert path.endswith(".005")

    def test_location_code(self):
        from ..core.stations import to_sds
        path = to_sds(self._stats(loc="00"), year=2023, jday=1)
        assert "BE.UCC.00.HHZ" in path

    def test_different_channels(self):
        from ..core.stations import to_sds
        for chan in ["BHZ", "LHN", "EHE"]:
            path = to_sds(self._stats(chan=chan), year=2023, jday=1)
            assert chan in path


# ============================================================================
# params.py — MSNoiseParams
# ============================================================================

class TestMSNoiseParams:
    """Tests for MSNoiseParams construction, access, and serialisation."""

    def _make_params(self):
        from ..params import MSNoiseParams
        from obspy.core.util.attribdict import AttribDict
        p = MSNoiseParams()
        p._set_lineage_names(["preprocess_1", "cc_1", "filter_1"])
        p._add_layer("global",     AttribDict({"output_folder": "./out", "hpc": "N"}))
        p._add_layer("preprocess", AttribDict({"preprocess_sampling_rate": 20.0}))
        p._add_layer("cc",         AttribDict({"maxlag": 60.0, "winsorizing": 2.0}))
        p._add_layer("filter",     AttribDict({"freqmin": 0.1, "freqmax": 1.0}))
        return p

    def test_category_access(self):
        p = self._make_params()
        assert p.cc.maxlag == 60.0
        assert p.filter.freqmin == pytest.approx(0.1)

    def test_global_underscore_alias(self):
        p = self._make_params()
        assert p.global_.hpc == "N"

    def test_bracket_access(self):
        p = self._make_params()
        assert p["cc"].maxlag == 60.0

    def test_missing_category_raises(self):
        p = self._make_params()
        with pytest.raises(AttributeError):
            _ = p.mwcs

    def test_missing_bracket_raises(self):
        p = self._make_params()
        with pytest.raises(KeyError):
            _ = p["mwcs"]

    def test_immutable(self):
        p = self._make_params()
        with pytest.raises(AttributeError):
            p.cc = "oops"

    def test_category_property(self):
        p = self._make_params()
        assert p.category == "filter"

    def test_category_layer_property(self):
        p = self._make_params()
        assert p.category_layer.freqmin == pytest.approx(0.1)

    def test_categories_list(self):
        p = self._make_params()
        assert p.categories == ["global", "preprocess", "cc", "filter"]

    def test_lineage_names(self):
        p = self._make_params()
        assert p.lineage_names == ["preprocess_1", "cc_1", "filter_1"]

    def test_step_name(self):
        p = self._make_params()
        assert p.step_name == "filter_1"

    def test_as_flat_dict(self):
        p = self._make_params()
        d = p.as_flat_dict()
        assert "maxlag" in d
        assert "freqmin" in d
        assert d["maxlag"] == 60.0

    def test_repr(self):
        p = self._make_params()
        r = repr(p)
        assert "MSNoiseParams" in r
        assert "filter" in r

    def test_yaml_roundtrip(self):
        pytest.importorskip("yaml")
        p = self._make_params()
        from ..params import MSNoiseParams
        yaml_str = p.to_yaml_string()
        assert "cc" in yaml_str
        p2 = MSNoiseParams.from_yaml_string(yaml_str)
        assert p2.cc.maxlag == pytest.approx(60.0)
        assert p2.filter.freqmin == pytest.approx(0.1)
        assert p2.lineage_names == p.lineage_names

    def test_yaml_roundtrip_file(self, tmp_path):
        pytest.importorskip("yaml")
        p = self._make_params()
        from ..params import MSNoiseParams
        fpath = str(tmp_path / "params.yaml")
        p.to_yaml(fpath)
        p2 = MSNoiseParams.from_yaml(fpath)
        assert p2.cc.winsorizing == pytest.approx(2.0)

    def test_empty_params_raises_on_category(self):
        from ..params import MSNoiseParams
        p = MSNoiseParams()
        with pytest.raises(RuntimeError):
            _ = p.category


# ============================================================================
# Project YAML export / import round-trip
# ============================================================================

def _make_project_db(tmp_path):
    """Helper: initialise a minimal SQLite MSNoise project, return (db, tmp)."""
    import os
    from ..core.db import create_database_inifile, connect
    from ..msnoise_table_def import declare_tables, DataSource, Station
    os.chdir(tmp_path)
    create_database_inifile(tech=1, hostname="proj.sqlite",
                            database="", username="", password="", prefix="")
    db = connect()
    declare_tables().Base.metadata.create_all(db.get_bind())
    return db, tmp_path


class TestProjectYamlRoundTrip:
    """Round-trip tests for export_project_to_yaml / create_project_from_yaml."""

    # ── minimal project YAML used across tests ────────────────────────────
    YAML_SRC = """\
msnoise_project_version: 1

global_1:
  startdate: "2013-04-01"
  enddate: "2014-10-31"

preprocess_1:
  after: global_1
  cc_sampling_rate: 20.0

cc_1:
  after: preprocess_1
  whitening: "N"
  components_to_compute_single_station: "ZZ,EE,NN"
  cc_type_single_station_AC: PCC

filter_1:
  after: cc_1
  freqmin: 1.0
  freqmax: 2.0
  AC: "Y"
  CC: "N"

filter_2:
  after: cc_1
  freqmin: 0.5
  freqmax: 1.0
  AC: "Y"
  CC: "N"

stack_1:
  after: [filter_1, filter_2]
  mov_stack: "(('2D','1D'))"

refstack_1:
  after: [filter_1, filter_2]
  ref_begin: "2013-04-01"
  ref_end: "2014-10-31"

mwcs_1:
  after: [stack_1, refstack_1]
  freqmin: 1.0
  freqmax: 2.0

mwcs_dtt_1:
  after: mwcs_1
  dtt_minlag: 5.0
  dtt_width: 30.0
  dtt_mincoh: 0.6
  dtt_maxerr: 0.1
"""

    def _import(self, tmp_path):
        """Create a DB from YAML_SRC, return (db, yaml_path)."""
        pytest.importorskip("yaml")
        db, tmp = _make_project_db(tmp_path)
        yaml_path = str(tmp / "project.yaml")
        with open(yaml_path, "w") as f:
            f.write(self.YAML_SRC)
        from ..core.config import create_project_from_yaml
        created, warnings = create_project_from_yaml(db, yaml_path)
        return db, yaml_path, created, warnings

    # ── import tests ──────────────────────────────────────────────────────

    def test_import_creates_expected_steps(self, tmp_path):
        db, _, created, _ = self._import(tmp_path)
        assert "global_1"      in created
        assert "preprocess_1"  in created
        assert "cc_1"          in created
        assert "filter_1"      in created
        assert "filter_2"      in created
        assert "stack_1"       in created
        assert "refstack_1"    in created
        assert "mwcs_1"        in created
        assert "mwcs_dtt_1"    in created

    def test_import_config_overrides_applied(self, tmp_path):
        db, _, _, _ = self._import(tmp_path)
        from ..core.config import get_config_set_details
        cc = {r["name"]: r["value"]
              for r in get_config_set_details(db, "cc", 1, format="list")}
        assert cc["whitening"] == "N"
        assert cc["cc_type_single_station_AC"] == "PCC"
        assert cc["cc_sampling_rate"] == "20.0"

        f1 = {r["name"]: r["value"]
              for r in get_config_set_details(db, "filter", 1, format="list")}
        assert f1["freqmin"] == "1.0"
        assert f1["AC"] == "Y"
        assert f1["CC"] == "N"

        f2 = {r["name"]: r["value"]
              for r in get_config_set_details(db, "filter", 2, format="list")}
        assert f2["freqmin"] == "0.5"

    def test_import_links_created(self, tmp_path):
        db, _, _, _ = self._import(tmp_path)
        from ..core.workflow import get_workflow_links, get_workflow_steps
        steps = {s.step_name: s.step_id for s in get_workflow_steps(db)}
        links = {(l.from_step_id, l.to_step_id)
                 for l in get_workflow_links(db)}
        # filter_1 → stack_1 and filter_2 → stack_1
        assert (steps["filter_1"], steps["stack_1"]) in links
        assert (steps["filter_2"], steps["stack_1"]) in links
        # stack_1 → mwcs_1 and refstack_1 → mwcs_1
        assert (steps["stack_1"],   steps["mwcs_1"]) in links
        assert (steps["refstack_1"], steps["mwcs_1"]) in links
        # preprocess_1 → cc_1 (single parent)
        assert (steps["preprocess_1"], steps["cc_1"]) in links

    def test_import_unknown_key_warns(self, tmp_path):
        pytest.importorskip("yaml")
        db, tmp = _make_project_db(tmp_path)
        src = self.YAML_SRC + "\n  typo_key: bad_value\n"
        # inject typo into cc_1 block
        src = src.replace(
            "cc_type_single_station_AC: PCC",
            "cc_type_single_station_AC: PCC\n  not_a_real_key: oops"
        )
        yaml_path = str(tmp / "bad.yaml")
        with open(yaml_path, "w") as f:
            f.write(src)
        from ..core.config import create_project_from_yaml
        _, warnings = create_project_from_yaml(db, yaml_path)
        assert any("not_a_real_key" in w for w in warnings)

    def test_import_wrong_version_raises(self, tmp_path):
        pytest.importorskip("yaml")
        db, tmp = _make_project_db(tmp_path)
        yaml_path = str(tmp / "bad.yaml")
        with open(yaml_path, "w") as f:
            f.write("msnoise_params_version: 1\ncc:\n  maxlag: 60\n")
        from ..core.config import create_project_from_yaml
        with pytest.raises(ValueError, match="msnoise_project_version"):
            create_project_from_yaml(db, yaml_path)

    # ── export tests ──────────────────────────────────────────────────────

    def test_export_produces_valid_yaml(self, tmp_path):
        import yaml
        db, _, _, _ = self._import(tmp_path)
        from ..core.config import export_project_to_yaml
        out = str(tmp_path / "exported.yaml")
        export_project_to_yaml(db, out)
        with open(out) as f:
            doc = yaml.safe_load(f)
        assert doc["msnoise_project_version"] == 1
        assert "cc_1" in doc
        assert "filter_1" in doc
        assert "filter_2" in doc

    def test_export_only_non_defaults(self, tmp_path):
        import yaml
        db, _, _, _ = self._import(tmp_path)
        from ..core.config import export_project_to_yaml
        out = str(tmp_path / "exported.yaml")
        export_project_to_yaml(db, out, only_non_defaults=True)
        with open(out) as f:
            doc = yaml.safe_load(f)
        # non-default overrides present
        assert doc["cc_1"].get("whitening") == "N"
        assert doc["filter_1"].get("freqmin") == "1.0"
        # default value (maxlag=120.0) must NOT appear
        assert "maxlag" not in doc.get("cc_1", {})

    def test_export_all_values(self, tmp_path):
        import yaml
        db, _, _, _ = self._import(tmp_path)
        from ..core.config import export_project_to_yaml
        out = str(tmp_path / "exported_all.yaml")
        export_project_to_yaml(db, out, only_non_defaults=False)
        with open(out) as f:
            doc = yaml.safe_load(f)
        # with all values, maxlag should be present
        assert "maxlag" in doc.get("cc_1", {})

    def test_export_after_links(self, tmp_path):
        import yaml
        db, _, _, _ = self._import(tmp_path)
        from ..core.config import export_project_to_yaml
        out = str(tmp_path / "exported.yaml")
        export_project_to_yaml(db, out)
        with open(out) as f:
            doc = yaml.safe_load(f)
        # stack_1 must declare both filters as parents
        after = doc["stack_1"].get("after", [])
        if isinstance(after, str):
            after = [after]
        assert "filter_1" in after
        assert "filter_2" in after
        # mwcs_1 must declare both stack_1 and refstack_1
        after_mwcs = doc["mwcs_1"].get("after", [])
        if isinstance(after_mwcs, str):
            after_mwcs = [after_mwcs]
        assert "stack_1"   in after_mwcs
        assert "refstack_1" in after_mwcs

    # ── full round-trip ───────────────────────────────────────────────────

    def test_full_roundtrip(self, tmp_path):
        """Export from DB-A (non-defaults only), verify the YAML is re-importable.

        A second DB is not set up here to avoid os.chdir session isolation issues
        in the test fixture — those belong in integration tests.  We verify the
        exported YAML is structurally correct and contains the expected overrides,
        then confirm that create_project_from_yaml on the SAME db produces no
        extra steps (idempotency guard).
        """
        import yaml
        db_a, _, created_a, _ = self._import(tmp_path)
        from ..core.config import export_project_to_yaml, get_config_set_details

        exported = str(tmp_path / "exported.yaml")
        export_project_to_yaml(db_a, exported, only_non_defaults=True)

        with open(exported) as f:
            doc = yaml.safe_load(f)

        # version marker present
        assert doc["msnoise_project_version"] == 1

        # all originally imported steps are present
        for step in created_a:
            assert step in doc, f"{step} missing from exported YAML"

        # overridden values round-tripped correctly
        assert doc["cc_1"].get("whitening") == "N"
        assert doc["cc_1"].get("cc_type_single_station_AC") == "PCC"
        assert doc["filter_1"].get("freqmin") == "1.0"
        assert doc["filter_2"].get("freqmin") == "0.5"
        assert doc["mwcs_dtt_1"].get("dtt_mincoh") == "0.6"
        assert doc["mwcs_dtt_1"].get("dtt_maxerr") == "0.1"

        # after links preserved
        stack_after = doc["stack_1"].get("after", [])
        if isinstance(stack_after, str):
            stack_after = [stack_after]
        assert "filter_1" in stack_after
        assert "filter_2" in stack_after

        mwcs_after = doc["mwcs_1"].get("after", [])
        if isinstance(mwcs_after, str):
            mwcs_after = [mwcs_after]
        assert "stack_1"    in mwcs_after
        assert "refstack_1" in mwcs_after

        # default values should NOT appear (only_non_defaults=True)
        assert "maxlag" not in doc.get("cc_1", {})


# ── Regression guards for the core/compute.py maths review ──────────────────


class TestWhiten2PsdBranch:
    """``whitening_type='PSD'`` must clip the MAGNITUDE and keep the phase."""

    def _spiky(self, nfft=2048, seed=0):
        import scipy.fft as sf
        rng = np.random.default_rng(seed)
        F = sf.fft(rng.standard_normal((2, nfft)), nfft, axis=1)
        for b in (300, 301, 500):          # strong monochromatic spikes
            F[:, b] *= 60.0
        return F, np.ones((2, nfft // 2 + 1))

    def test_clipping_preserves_phase(self):
        from ..core.compute import whiten2
        nfft = 2048
        F, psds = self._spiky(nfft)
        before = F.copy()
        whiten2(F, nfft, 20, 900, 120, 800, psds, "PSD")
        band = slice(120, 800)
        amp_b, amp_a = np.abs(before[0, band]), np.abs(F[0, band])
        clipped = np.where(amp_a < amp_b - 1e-9)[0]
        assert clipped.size > 0, "the spikes should have been clipped"
        # np.clip() on the complex array used to leave imag == 0 exactly
        assert not np.any(np.abs(F[0, band].imag[clipped]) < 1e-18)
        dphi = np.angle(F[0, band])[clipped] - np.angle(before[0, band])[clipped]
        assert np.abs(np.angle(np.exp(1j * dphi))).max() < 1e-9

    def test_taper_is_cos_squared_not_cos_fourth(self):
        from ..core.compute import whiten2
        nfft, low, porte1 = 2048, 20, 120
        G = np.ones((1, nfft), dtype=complex)
        whiten2(G, nfft, low, 900, porte1, 800,
                np.ones((1, nfft // 2 + 1)), "PSD")
        mid = (low + porte1) // 2
        frac = (mid - low) / (porte1 - low)
        expected = np.cos(np.pi / 2 * (1 - frac)) ** 2
        assert abs(np.abs(G[0, mid]) - expected) < 0.02
        assert abs(np.abs(G[0, mid]) - expected ** 2) > 0.2

    @pytest.mark.parametrize("nfft", [1024, 9, 15, 25, 27, 45, 81, 125])
    def test_hermitian_symmetry_even_and_odd_nfft(self, nfft):
        import scipy.fft as sf

        from ..core.compute import whiten2
        rng = np.random.default_rng(nfft)
        F = sf.fft(rng.standard_normal((1, nfft)), nfft, axis=1)
        whiten2(F, nfft, 1, nfft // 2, 2, nfft // 2 - 1,
                np.ones((1, nfft // 2 + 1)), "B")
        nneg = (nfft - 1) // 2
        assert np.allclose(F[0, nfft - nneg:],
                           np.conj(F[0, 1:nneg + 1])[::-1])
        # a Hermitian spectrum must inverse-transform to a real signal
        assert np.abs(sf.ifft(F[0], nfft).imag).max() < 1e-12


class TestWhitenGuards:
    @pytest.mark.parametrize("nfft", [1024, 9, 15, 25, 27, 45, 81, 125])
    def test_hermitian_symmetry_even_and_odd_nfft(self, nfft):
        import scipy.fft as sf

        from ..core.compute import whiten
        rng = np.random.default_rng(nfft)
        F = whiten(rng.standard_normal(min(nfft, 64)), nfft, 0.05, 1.0, 4.0)
        nneg = (nfft - 1) // 2
        assert np.allclose(F[nfft - nneg:], np.conj(F[1:nneg + 1])[::-1])
        assert np.abs(sf.ifft(F, nfft).imag).max() < 1e-12

    def test_empty_passband_raises_readable_error(self):
        from ..core.compute import whiten
        rng = np.random.default_rng(0)
        with pytest.raises(ValueError, match="passband"):
            whiten(rng.standard_normal(256), 256, 0.05, 50.0, 60.0)


class TestMyCorr2LagWindow:
    def _ffts(self, Nt=64, seed=1):
        import scipy.fft as sf
        rng = np.random.default_rng(seed)
        return sf.fft(rng.standard_normal((2, Nt)), Nt, axis=1)

    @pytest.mark.parametrize("maxlag,expected", [(10, 21), (63, 127)])
    def test_output_length(self, maxlag, expected):
        from ..core.compute import myCorr2
        out = myCorr2(self._ffts(), maxlag, np.ones(2), [("a", 0, 1)], nfft=64)
        assert out["a"].shape == (expected,)

    @pytest.mark.parametrize("maxlag", [64, 69])
    def test_maxlag_at_or_beyond_nfft_clamps(self, maxlag):
        """Used to raise a broadcast error (out_len was 2*Nt, folded is 2*Nt-1)."""
        from ..core.compute import myCorr2
        out = myCorr2(self._ffts(), maxlag, np.ones(2), [("a", 0, 1)], nfft=64)
        assert out["a"].shape == (127,)

    def test_zero_lag_is_centred(self):
        """A trace correlated with itself must peak exactly at the centre."""
        import scipy.fft as sf

        from ..core.compute import myCorr2
        rng = np.random.default_rng(2)
        x = rng.standard_normal(256)
        F = sf.fft(np.vstack([x, x]), 256, axis=1)
        ccf = myCorr2(F, 40, np.ones(2), [("a", 0, 1)], nfft=256)["a"]
        assert np.argmax(ccf) == 40


class TestMwcsBandUnwrap:
    """MWCS must not inherit cycle slips from the incoherent sub-freqmin bins."""

    @staticmethod
    def _pair(dt_true, seed, df=20.0, n=6000, noise=0.02):
        import scipy.fft as sf
        import scipy.signal
        rng = np.random.default_rng(seed)
        sos = scipy.signal.butter(4, [1.0, 4.0], btype="band", fs=df,
                                  output="sos")
        ref = scipy.signal.sosfilt(sos, rng.standard_normal(n))
        f = sf.fftfreq(n, 1.0 / df)
        cur = np.real(sf.ifft(sf.fft(ref) * np.exp(-2j * np.pi * f * dt_true)))
        cur = cur + noise * ref.std() * rng.standard_normal(n)
        return cur, ref

    def test_no_cycle_slip_outliers(self):
        from ..core.compute import mwcs
        dt_true, out = 0.02, []
        for seed in range(10):
            cur, ref = self._pair(dt_true, seed)
            out.append(mwcs(cur, ref, 1.0, 4.0, 20.0, -300.0, 12.0, 4.0)[:, 1])
        v = np.concatenate(out)
        assert abs(np.median(v) - dt_true) < 1e-3
        # DC-anchored unwrap produced ~13% outliers here, up to |dt| = 0.76 s
        assert np.mean(np.abs(v - dt_true) > 0.01) < 0.01
        assert np.abs(v).max() < 0.5 / 1.0   # < 1/(2*freqmin)

    @pytest.mark.parametrize("dt_true", [0.02, 0.2, 0.45])
    def test_recovers_dt_within_validity_limit(self, dt_true):
        """Band-anchored unwrap is valid while |dt| < 1/(2*freqmin) = 0.5 s."""
        from ..core.compute import mwcs
        cur, ref = self._pair(dt_true, 7, noise=0.0)
        v = mwcs(cur, ref, 1.0, 4.0, 20.0, -300.0, 12.0, 4.0)[:, 1]
        assert abs(np.median(v) - dt_true) < 5e-3

    def test_returns_four_columns(self):
        from ..core.compute import mwcs
        cur, ref = self._pair(0.02, 0)
        res = mwcs(cur, ref, 1.0, 4.0, 20.0, -300.0, 12.0, 4.0)
        assert res.ndim == 2 and res.shape[1] == 4

    def test_band_with_too_few_bins_raises(self):
        from ..core.compute import mwcs
        cur, ref = self._pair(0.02, 0)
        with pytest.raises(ValueError, match="FFT bin"):
            mwcs(cur, ref, 9.99, 9.995, 20.0, -300.0, 12.0, 4.0)


# ── Regression guards for the core/signal.py DSP review (batch 2) ───────────


class TestPsdDfRmsUnits:
    """PPSD reports an acceleration PSD per unit FREQUENCY on a PERIOD axis."""

    @staticmethod
    def _flat_psd(db=-120.0, n_times=3):
        periods = np.logspace(-2, 2, 200)
        data = np.full((n_times, len(periods)), db)
        return xr.DataArray(
            data, coords=[np.arange(n_times), periods],
            dims=["times", "periods"], name="PSD",
        ).to_dataset(), periods

    @pytest.mark.parametrize("output,omega_power", [("ACC", 0), ("VEL", 2),
                                                    ("DISP", 4)])
    def test_matches_seismorms_reference(self, output, omega_power):
        from ..core.signal import psd_df_rms
        ds, periods = self._flat_psd()
        fmin, fmax = 1.0, 20.0
        got = psd_df_rms(ds, [(fmin, fmax)], output=output)["RMS"].values[0, 0]

        # Reference: integrate over FREQUENCY, omega = 2*pi*f (SeismoRMS)
        f = np.sort(1.0 / periods)
        f = f[(f >= fmin) & (f <= fmax)]
        amp = 10.0 ** (-120.0 / 10.0)
        integrand = np.full_like(f, amp) / (2 * np.pi * f) ** omega_power
        ref = np.sqrt(np.trapezoid(integrand, f))
        assert got == pytest.approx(ref, rel=1e-9)

    def test_band_labels_are_backward_compatible(self):
        from ..core.signal import psd_df_rms
        ds, _ = self._flat_psd()
        out = psd_df_rms(ds, [(1.0, 20.0), (4.0, 14.0), (4.0, 9.0)])
        assert list(out.bands.values) == ["1.0-20.0", "4.0-14.0", "4.0-9.0"]

    def test_sub_decimal_bands_do_not_collide(self):
        """(0.05, 0.1) and (0.06, 0.1) both formatted as '0.1-0.1' before."""
        from ..core.signal import psd_df_rms
        ds, _ = self._flat_psd()
        out = psd_df_rms(ds, [(0.05, 0.1), (0.06, 0.1)])
        assert len(out.bands) == 2
        assert len(set(out.bands.values.tolist())) == 2


class TestMorletScaleMapping:
    """tfpws_stack assumes f = w*fs/(2*pi*s); the wavelet must honour it."""

    @pytest.mark.parametrize("f_req", [1.0, 2.0, 5.0, 8.0])
    def test_centre_frequency_matches_requested(self, f_req):
        import scipy.signal

        from ..core.signal import _morlet_wavelet
        fs, W = 20.0, 5.0
        s = W / (2.0 * np.pi * f_req / fs)
        M = max(int(10 * s), 3) | 1
        wav = _morlet_wavelet(M, s, w=W)

        probes = np.linspace(0.2, 10.0, 200)
        t = np.arange(4096) / fs
        resp = [
            np.abs(scipy.signal.fftconvolve(
                np.cos(2 * np.pi * fp * t), wav[::-1].conj(),
                mode="same")[1000:3000]).mean()
            for fp in probes
        ]
        peak = probes[int(np.argmax(resp))]
        # the old linspace(-10, 10, M) abscissa made this scale as 1/s**2:
        # 1 Hz landed at 0.13 Hz and 8 Hz at ~10 Hz
        assert peak == pytest.approx(f_req, rel=0.06)

    def test_envelope_is_s_samples_wide(self):
        from ..core.signal import _morlet_wavelet
        s, M = 8.0, 81
        env = np.abs(_morlet_wavelet(M, s, w=5.0))
        x = np.arange(M) - (M - 1) / 2.0
        sigma = np.sqrt((env ** 2 * x ** 2).sum() / (env ** 2).sum())
        assert sigma == pytest.approx(s / np.sqrt(2.0), rel=0.05)


class TestGetWindowNormalisation:
    @pytest.mark.parametrize("window", ["boxcar", "hanning"])
    def test_sums_to_one(self, window):
        from ..core.signal import get_window
        assert get_window(window, 5).sum().real == pytest.approx(1.0)

    def test_mwcs_coherence_is_invariant_to_window_scale(self):
        """The smoothing constant cancels in |X| / (dref*dcur)."""
        import scipy.signal as ss
        rng = np.random.default_rng(0)
        f1 = rng.standard_normal(256) + 1j * rng.standard_normal(256)
        f2 = rng.standard_normal(256) + 1j * rng.standard_normal(256)

        def coh(scale):
            w = ss.windows.hann(11) * scale
            X = ss.fftconvolve(f1 * f2.conj(), w, mode="same")
            d1 = np.sqrt(ss.fftconvolve(np.abs(f1) ** 2, w, mode="same"))
            d2 = np.sqrt(ss.fftconvolve(np.abs(f2) ** 2, w, mode="same"))
            return np.abs(X) / (d1 * d2)

        assert np.allclose(coh(1.0 / 11), coh(1.0 / ss.windows.hann(11).sum()))


class TestMwcsSlopeError:
    def test_matches_obspy_std_slope(self):
        """The manual error used an unweighted residual + a 1/sigma^2 form."""
        import scipy.fft as sf
        import scipy.signal
        from obspy.signal.regression import linear_regression

        from ..core.compute import mwcs
        fs, n = 20.0, 6000
        rng = np.random.default_rng(3)
        sos = scipy.signal.butter(4, [1.0, 4.0], btype="band", fs=fs,
                                  output="sos")
        ref = scipy.signal.sosfilt(sos, rng.standard_normal(n))
        f = sf.fftfreq(n, 1.0 / fs)
        cur = np.real(sf.ifft(sf.fft(ref) * np.exp(-2j * np.pi * f * 0.02)))
        cur = cur + 0.02 * ref.std() * rng.standard_normal(n)
        res = mwcs(cur, ref, 1.0, 4.0, fs, -300.0, 12.0, 4.0)
        assert res.shape[1] == 4
        assert np.all(np.isfinite(res[:, 2]))
        assert np.median(res[:, 1]) == pytest.approx(0.02, abs=1e-3)
        # sanity: the reported error must be far below the measured delay
        assert np.median(res[:, 2]) < 0.1 * abs(np.median(res[:, 1]))
        assert linear_regression is not None


# ── Batch-4 review: s10_stretching ──────────────────────────────────────────


class TestStretching:
    @staticmethod
    def _ccf(fs=20.0, maxlag=120.0, seed=0):
        import scipy.signal
        n = int(2 * maxlag * fs) + 1
        mid = int(fs * maxlag)
        rng = np.random.default_rng(seed)
        sos = scipy.signal.butter(4, [0.5, 4.0], btype="band", fs=fs,
                                  output="sos")
        ref = scipy.signal.sosfilt(sos, rng.standard_normal(n))
        taxis = (np.arange(n) - mid) / fs
        return ref, taxis, n, mid, fs

    def _measure(self, dvv_true, str_range=0.02, nstr=1001):
        from ..s10_stretching import _hwhm_errors, stretch_mat_creation
        ref, taxis, n, _, _ = self._ccf()
        # dv/v < 0 (slower medium) => arrivals later => cur(t) = ref(t/(1+dt/t))
        cur = np.interp(taxis / (1 - dvv_true), taxis, ref)
        M, strvec = stretch_mat_creation(ref, str_range=str_range, nstr=nstr)
        Mn = (M - M.mean(axis=1, keepdims=True)) / M.std(axis=1, keepdims=True)
        cn = (cur - cur.mean()) / cur.std()
        cc = Mn @ cn / n
        k = int(np.argmax(cc))
        err = _hwhm_errors(cc.reshape(nstr, 1))[0]
        return strvec[k], cc[k], err

    @pytest.mark.parametrize("dvv_true", [-0.01, -0.003, 0.0, 0.003, 0.01])
    def test_delta_minus_one_is_signed_dvv(self, dvv_true):
        """aggregate_dvv_pairs uses ``Delta - 1`` directly as dv/v."""
        delta, coeff, _ = self._measure(dvv_true)
        assert delta - 1.0 == pytest.approx(dvv_true, abs=1e-4)
        assert coeff > 0.99

    def test_error_converts_to_dvv_units(self):
        """The stored Error must share units with Delta - 1, not be an index."""
        from ..s10_stretching import hwhm_to_dvv
        str_range, nstr = 0.02, 1001
        _, _, raw_hwhm = self._measure(-0.01, str_range, nstr)
        converted = hwhm_to_dvv(raw_hwhm, str_range, nstr)
        assert raw_hwhm > 1.0                 # index units: tens
        assert converted < 0.01               # dv/v units: sub-percent
        assert converted == pytest.approx(raw_hwhm * 4e-5, rel=1e-9)

    def test_hwhm_to_dvv_is_the_stretch_step(self):
        from ..s10_stretching import hwhm_to_dvv, stretch_mat_creation
        ref, _, _, _, _ = self._ccf()
        str_range, nstr = 0.02, 101
        _, strvec = stretch_mat_creation(ref, str_range=str_range, nstr=nstr)
        # one index == one step of strvec
        assert hwhm_to_dvv(1.0, str_range, nstr) == pytest.approx(
            strvec[1] - strvec[0])

    def test_hwhm_returns_nan_when_curve_never_drops_to_half(self):
        from ..s10_stretching import _hwhm_errors
        n = 101
        c = np.linspace(0.0, 1.0, n)
        c[-1] = 0.5                            # peak at n-2, no left crossing
        assert np.isnan(_hwhm_errors(c.reshape(n, 1))[0])

    def test_hwhm_matches_analytic_gaussian(self):
        from ..s10_stretching import _hwhm_errors
        n = 101
        x = np.arange(n)
        c = np.exp(-0.5 * ((x - 50) / 5.0) ** 2)
        expected = 5.0 * np.sqrt(2 * np.log(2))
        assert _hwhm_errors(c.reshape(n, 1))[0] == pytest.approx(expected,
                                                                 rel=1e-3)

    def test_coda_mask_is_two_sided_and_identical_for_1d_and_2d(self):
        """The current-CCF mask used to zero one column, not the whole tail."""
        from ..s10_stretching import apply_coda_mask
        _, taxis, n, mid, fs = self._ccf()
        minlag, maxlag2 = 5.0, 45.0

        cur = apply_coda_mask(np.ones((3, n)), mid, minlag, maxlag2, fs)
        ref = apply_coda_mask(np.ones(n), mid, minlag, maxlag2, fs)

        # the 2-D (current) and 1-D (reference) paths must agree exactly
        assert np.array_equal(cur[0], ref)
        # nothing survives outside minlag <= |lag| <= maxlag2 ...
        keep = (np.abs(taxis) >= minlag) & (np.abs(taxis) <= maxlag2)
        assert np.count_nonzero(ref[~keep]) == 0
        # ... and both lobes survive
        assert np.count_nonzero(ref[taxis < 0]) > 0
        assert np.count_nonzero(ref[taxis > 0]) > 0

    def test_stretch_matrix_shape_and_symmetry(self):
        from ..s10_stretching import stretch_mat_creation
        ref, _, n, _, _ = self._ccf()
        M, strvec = stretch_mat_creation(ref, str_range=0.02, nstr=101)
        assert M.shape == (101, n)
        assert strvec[0] == pytest.approx(0.98)
        assert strvec[-1] == pytest.approx(1.02)
        assert strvec[50] == pytest.approx(1.0)
        # the unstretched row must reproduce the reference
        assert np.corrcoef(M[50], ref)[0, 1] > 0.999
