"""Regression tests for the self-contained wavelet implementation in
``msnoise.core.signal``.

Reference arrays in ``fixtures/pycwt_reference.npz`` were generated with
pycwt 0.4.0b1 (the version previously required by MSNoise) using the
signals and parameters embedded in that file.  The tests below verify
that the internal reimplementation produces numerically identical results
(to within floating-point rounding, <= 1e-10) for every code path that
was previously delegated to pycwt:

- :class:`_Morlet` - default and non-default ``f0`` (production path)
- :class:`_Paul`   - order 4 (T&C 1998 Table 1 row 2)
- :class:`_DOG`    - order 2 = Mexican hat (Table 1 row 3)
- :func:`_cwt`     - auto-scale mode AND custom-frequency mode
- Full WCT pipeline: :func:`get_wavelet_type` -> :func:`prepare_ref_wct`
  -> :func:`apply_wct` (shape / finite / monotone-freq checks)
- :func:`stack` with ``stack_method="tfpws"`` (Schimmel & Gallart 2007)
"""
import warnings
from pathlib import Path

import numpy as np
import pytest

FIXTURE = Path(__file__).parent / "fixtures" / "pycwt_reference.npz"


def _load():
    if not FIXTURE.exists():
        pytest.skip("pycwt_reference.npz fixture not found")
    return np.load(FIXTURE)


class TestMorletWavelet:
    def test_psi_ft_peak_at_f0(self):
        from ..core.signal import _Morlet
        m = _Morlet(6)
        f = np.linspace(0, 15, 500)
        vals = np.abs(m.psi_ft(f))
        assert f[np.argmax(vals)] == pytest.approx(6.0, abs=0.1)

    def test_flambda_positive(self):
        from ..core.signal import _Morlet
        assert _Morlet(6).flambda() > 0

    def test_coi_positive(self):
        from ..core.signal import _Morlet
        assert _Morlet(6).coi() > 0

    def test_f0_parameter(self):
        from ..core.signal import _Morlet
        m4 = _Morlet(4)
        f = np.linspace(0, 10, 200)
        assert f[np.argmax(np.abs(m4.psi_ft(f)))] == pytest.approx(4.0, abs=0.1)


class TestPaulWavelet:
    def test_psi_ft_zero_for_negative_f(self):
        from ..core.signal import _Paul
        f = np.linspace(-5, 0, 50)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            assert np.allclose(_Paul(4).psi_ft(f), 0.0)

    def test_psi_ft_norm_matches_pycwt(self):
        """Paul psi_ft constant matches pycwt normalization (off-by-one fix)."""
        from ..core.signal import _Paul, _cwt
        ref = _load()
        sig = ref["sig_512"]
        dt, dj, s0 = float(ref["dt"]), float(ref["dj"]), float(ref["s0"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            W, *_ = _cwt(sig, dt, dj, s0, -1, _Paul(4))
        W_ref = ref["W_p4"]
        assert W.shape == W_ref.shape
        assert np.abs(np.abs(W) - np.abs(W_ref)).max() < 1e-10


class TestDOGWavelet:
    def test_psi_ft_complex(self):
        from ..core.signal import _DOG
        f = np.linspace(0, 5, 50)
        result = _DOG(2).psi_ft(f)
        assert np.iscomplexobj(result)

    def test_mexican_hat_is_dog2(self):
        from ..core.signal import _DOG, _MexicanHat
        f = np.linspace(0, 5, 50)
        assert np.allclose(_MexicanHat().psi_ft(f), _DOG(2).psi_ft(f))


class TestCWT:
    """Full array regression: new _cwt must match pycwt output to 1e-10."""

    @pytest.fixture(autouse=True)
    def ref(self):
        self._ref = _load()

    def _run(self, sig_key, wavelet, freqs_key=None,
             W_key=None, sj_key=None, freqs_out_key=None, coi_key=None):
        from ..core.signal import _cwt
        ref = self._ref
        dt, dj, s0 = float(ref["dt"]), float(ref["dj"]), float(ref["s0"])
        sig   = ref[sig_key]
        freqs = ref[freqs_key] if freqs_key else None
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            W, sj, fr, coi, _, _ = _cwt(sig, dt, dj, s0, -1, wavelet, freqs=freqs)
        if W_key:
            assert W.shape == ref[W_key].shape
            err = np.abs(np.abs(W) - np.abs(ref[W_key])).max()
            assert err < 1e-10, f"|W| max error = {err:.2e}"
        if sj_key:
            assert np.allclose(sj, ref[sj_key], atol=0, rtol=1e-12)
        if freqs_out_key:
            assert np.allclose(fr, ref[freqs_out_key], atol=0, rtol=1e-12)
        if coi_key:
            assert np.allclose(coi, ref[coi_key], atol=0, rtol=1e-12)
        return W, sj, fr, coi

    def test_morlet6_auto_scales_amplitude(self):
        from ..core.signal import _Morlet
        self._run("sig_512", _Morlet(6),
                  W_key="W_m6", sj_key="sj_m6",
                  freqs_out_key="freqs_m6", coi_key="coi_m6")

    def test_morlet6_custom_freqs_amplitude(self):
        from ..core.signal import _Morlet
        self._run("sig_500", _Morlet(6), freqs_key="freqlim",
                  W_key="W_cf", sj_key="sj_cf",
                  freqs_out_key="freqs_cf", coi_key="coi_cf")

    def test_paul4_amplitude(self):
        from ..core.signal import _Paul
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self._run("sig_512", _Paul(4),
                      W_key="W_p4", sj_key="sj_p4",
                      freqs_out_key="freqs_p4", coi_key="coi_p4")

    def test_dog2_amplitude(self):
        from ..core.signal import _DOG
        self._run("sig_512", _DOG(2),
                  W_key="W_d2", sj_key="sj_d2",
                  freqs_out_key="freqs_d2", coi_key="coi_d2")

    def test_chirp_signal_amplitude(self):
        """Frequency-varying signal: validates time-resolution of the CWT."""
        from ..core.signal import _Morlet
        self._run("sig_chirp", _Morlet(6),
                  W_key="W_ch", sj_key="sj_ch",
                  freqs_out_key="freqs_ch", coi_key="coi_ch")

    def test_nonpow2_length_shape(self):
        """Non-power-of-two input must produce output trimmed to input length."""
        from ..core.signal import _cwt, _Morlet
        ref = self._ref
        sig = ref["sig_500"]
        W, *_ = _cwt(sig, float(ref["dt"]), float(ref["dj"]),
                     float(ref["s0"]), -1, _Morlet(6))
        assert W.shape[1] == len(sig)

    def test_coi_length_matches_signal(self):
        from ..core.signal import _cwt, _Morlet
        ref = self._ref
        sig = ref["sig_512"]
        _, _, _, coi, _, _ = _cwt(sig, float(ref["dt"]), float(ref["dj"]),
                                  float(ref["s0"]), -1, _Morlet(6))
        assert len(coi) == len(sig)

    def test_freqs_descending_auto(self):
        """Auto-scale frequencies must be strictly decreasing (coarse to fine)."""
        from ..core.signal import _cwt, _Morlet
        ref = self._ref
        _, _, freqs, *_ = _cwt(ref["sig_512"], float(ref["dt"]),
                               float(ref["dj"]), float(ref["s0"]),
                               -1, _Morlet(6))
        assert np.all(np.diff(freqs) < 0)


class TestGetWaveletType:
    def test_morlet_default_param(self):
        from ..core.signal import get_wavelet_type, _Morlet
        w = get_wavelet_type(("Morlet",))
        assert isinstance(w, _Morlet)
        assert w.f0 == 6

    def test_morlet_custom_param(self):
        from ..core.signal import get_wavelet_type, _Morlet
        w = get_wavelet_type(("Morlet", 4))
        assert isinstance(w, _Morlet) and w.f0 == 4.0

    def test_paul(self):
        from ..core.signal import get_wavelet_type, _Paul
        w = get_wavelet_type(("Paul", 4))
        assert isinstance(w, _Paul) and w.m == 4

    def test_dog(self):
        from ..core.signal import get_wavelet_type, _DOG
        w = get_wavelet_type(("DOG", 2))
        assert isinstance(w, _DOG) and w.m == 2

    def test_mexican_hat(self):
        from ..core.signal import get_wavelet_type, _MexicanHat
        w = get_wavelet_type(("MexicanHat",))
        assert isinstance(w, _MexicanHat)

    def test_unknown_raises(self):
        from ..core.signal import get_wavelet_type
        with pytest.raises(ValueError, match="Unknown wavelet type"):
            get_wavelet_type(("Bogus",))


class TestWCTPipeline:
    @pytest.fixture
    def signals(self):
        rng = np.random.default_rng(42)
        return rng.standard_normal(512), rng.standard_normal(512)

    def test_output_shape(self, signals):
        from ..core.signal import get_wavelet_type, prepare_ref_wct, apply_wct
        ref_sig, cur_sig = signals
        mother   = get_wavelet_type(("Morlet", 6))
        ref_data = prepare_ref_wct(ref_sig, fs=50., freqmin=1., freqmax=20.,
                                   nptsfreq=30, mother=mother)
        WXamp, WXspec, WXangle, Wcoh, WXdt, freqs, coi = apply_wct(ref_data, cur_sig)
        assert WXamp.shape == (30, 512)
        assert len(freqs)  == 30
        assert len(coi)    == 512

    def test_output_finite(self, signals):
        from ..core.signal import get_wavelet_type, prepare_ref_wct, apply_wct
        ref_sig, cur_sig = signals
        mother   = get_wavelet_type(("Morlet", 6))
        ref_data = prepare_ref_wct(ref_sig, fs=50., freqmin=1., freqmax=20.,
                                   nptsfreq=20, mother=mother)
        WXamp, _, _, Wcoh, *_ = apply_wct(ref_data, cur_sig)
        assert np.all(np.isfinite(WXamp))
        assert np.all(np.isfinite(Wcoh))

    def test_coherence_bounded(self, signals):
        from ..core.signal import get_wavelet_type, prepare_ref_wct, apply_wct
        ref_sig, cur_sig = signals
        mother   = get_wavelet_type(("Morlet", 6))
        ref_data = prepare_ref_wct(ref_sig, fs=50., freqmin=1., freqmax=20.,
                                   nptsfreq=20, mother=mother)
        *_, Wcoh, _, _, _ = apply_wct(ref_data, cur_sig)
        assert np.abs(Wcoh).max() <= 1.0 + 1e-6

    def test_freqs_monotone_decreasing(self, signals):
        from ..core.signal import get_wavelet_type, prepare_ref_wct, apply_wct
        ref_sig, cur_sig = signals
        mother   = get_wavelet_type(("Morlet", 6))
        ref_data = prepare_ref_wct(ref_sig, fs=50., freqmin=1., freqmax=20.,
                                   nptsfreq=20, mother=mother)
        *_, freqs, coi = apply_wct(ref_data, cur_sig)
        assert np.all(np.diff(freqs) < 0)

    def test_default_mother_is_morlet6(self):
        from ..core.signal import _Morlet, prepare_ref_wct
        rng      = np.random.default_rng(7)
        ref_data = prepare_ref_wct(rng.standard_normal(256), fs=50.,
                                   freqmin=1., freqmax=20., nptsfreq=10,
                                   mother=None)
        assert ref_data is not None
        assert isinstance(ref_data[8], _Morlet)
        assert ref_data[8].f0 == 6

    @pytest.mark.parametrize("wtype", [("Morlet", 6), ("Morlet", 4),
                                       ("Paul", 4), ("DOG", 2)])
    def test_mother_is_carried_to_apply_wct(self, wtype):
        """apply_wct must use the SAME mother as prepare_ref_wct.

        Regression guard: apply_wct previously hardcoded ``_Morlet(6)`` for
        the current-day CWT.  Since scales are ``1/(mother.flambda()*f)``, a
        mismatch puts cwt_ref and cwt_cur on different scale grids, so a
        trace correlated with itself no longer gives coherence 1 / zero lag.
        """
        from ..core.signal import apply_wct, get_wavelet_type, prepare_ref_wct
        rng    = np.random.default_rng(11)
        sig    = rng.standard_normal(512)
        mother = get_wavelet_type(wtype)
        ref_data = prepare_ref_wct(sig, fs=50., freqmin=1., freqmax=20.,
                                   nptsfreq=20, mother=mother)
        assert ref_data[8] is mother
        _, _, _, Wcoh, WXdt, _, _ = apply_wct(ref_data, sig)
        # self-coherence is identically 1 and the lag identically 0 only if
        # both CWTs used the same mother wavelet
        assert np.nanmin(np.abs(Wcoh)) == pytest.approx(1.0, abs=1e-6)
        assert np.nanmax(np.abs(WXdt)) == pytest.approx(0.0, abs=1e-9)


class TestTfPWSStack:
    @pytest.fixture
    def traces(self):
        return np.random.default_rng(99).standard_normal((25, 1000))

    def test_output_shape(self, traces):
        from ..core.signal import tfpws_stack
        assert tfpws_stack(traces, fs=50., freqmin=3., freqmax=12.).shape == (1000,)

    def test_output_real(self, traces):
        from ..core.signal import tfpws_stack
        assert not np.iscomplexobj(tfpws_stack(traces, fs=50., freqmin=3., freqmax=12.))

    def test_weight_suppresses_incoherent(self, traces):
        from ..core.signal import tfpws_stack
        lin  = traces.mean(axis=0)
        tfpw = tfpws_stack(traces, fs=50., freqmin=3., freqmax=12., power=2)
        assert np.abs(tfpw).max() <= np.abs(lin).max() + 1e-9

    def test_power_zero_equals_linear(self, traces):
        from ..core.signal import tfpws_stack
        lin  = traces.mean(axis=0)
        tfpw = tfpws_stack(traces, fs=50., freqmin=3., freqmax=12., power=0)
        assert np.allclose(tfpw, lin, atol=1e-12)

    def test_single_trace(self):
        from ..core.signal import tfpws_stack
        data   = np.random.default_rng(1).standard_normal((1, 500))
        result = tfpws_stack(data, fs=50., freqmin=3., freqmax=12.)
        assert np.allclose(result, data[0], atol=1e-12)

    def test_via_stack_dispatch(self, traces):
        from ..core.signal import stack, tfpws_stack
        direct   = tfpws_stack(traces, fs=50., freqmin=3., freqmax=12., power=2)
        via_disp = stack(traces, stack_method="tfpws",
                         goal_sampling_rate=50., freqmin=3., freqmax=12.,
                         pws_power=2, tfpws_nscales=20)
        assert np.allclose(direct, via_disp, atol=1e-12)

    def test_unknown_stack_method_falls_back_to_linear(self, traces):
        from ..core.signal import stack
        lin   = traces[~np.isnan(traces).any(axis=1)].mean(axis=0)
        bogus = stack(traces, stack_method="bogus_method", goal_sampling_rate=50.)
        assert np.allclose(bogus, lin, atol=1e-12)
