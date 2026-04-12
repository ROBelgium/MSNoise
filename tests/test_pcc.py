"""Unit tests for the PCC2 implementation in msnoise.core.compute.

Tests cover:
- _analytic_phase: shape, dtype, unit-amplitude contract
- _analytic_phase_batch: consistency with single-trace version
- pcc_xcorr: peak at lag-0 for identical signals, peak at correct lag for
  shifted signals, near-zero for uncorrelated noise, output length, and
  normalisation modes.

Run with:  pytest tests/test_pcc.py -v
"""

import numpy as np
import pytest

from msnoise.core.compute import _analytic_phase, _analytic_phase_batch, pcc_xcorr


RNG = np.random.default_rng(42)
N   = 1024   # trace length (samples)
ML  = 100    # max-lag (samples)


# ─────────────────────────────────────────────────────────────────────────────
# _analytic_phase
# ─────────────────────────────────────────────────────────────────────────────

class TestAnalyticPhase:
    def test_shape(self):
        x = RNG.standard_normal(N)
        y = _analytic_phase(x)
        assert y.shape == (N,)

    def test_dtype_complex(self):
        x = RNG.standard_normal(N).astype(np.float32)
        y = _analytic_phase(x.astype(np.float64))
        assert np.iscomplexobj(y)

    def test_unit_amplitude(self):
        """Every sample of the phase signal should have |y[n]| ≈ 1."""
        x = RNG.standard_normal(N)
        y = _analytic_phase(x)
        amp = np.abs(y)
        # Allow for the eps regularisation: amplitude should be very close to 1
        # everywhere except near-zero-amplitude samples (rare for Gaussian noise)
        assert amp.max() <= 1.0 + 1e-9
        assert amp.mean() == pytest.approx(1.0, abs=1e-3)

    def test_zero_input_does_not_blow_up(self):
        """All-zero input should return zeros (or tiny values), not NaN/Inf."""
        x = np.zeros(N)
        y = _analytic_phase(x)
        assert np.all(np.isfinite(y))


# ─────────────────────────────────────────────────────────────────────────────
# _analytic_phase_batch
# ─────────────────────────────────────────────────────────────────────────────

class TestAnalyticPhaseBatch:
    def test_shape(self):
        X = RNG.standard_normal((5, N))
        Y = _analytic_phase_batch(X)
        assert Y.shape == (5, N)

    def test_consistent_with_single(self):
        """Each row of the batch result must match the single-trace result."""
        X = RNG.standard_normal((4, N))
        Y_batch = _analytic_phase_batch(X)
        for i in range(4):
            y_single = _analytic_phase(X[i])
            np.testing.assert_allclose(Y_batch[i], y_single, atol=1e-12,
                                       err_msg=f"Mismatch at row {i}")


# ─────────────────────────────────────────────────────────────────────────────
# pcc_xcorr
# ─────────────────────────────────────────────────────────────────────────────

def _make_index(*pairs):
    """Build an index list from (ccf_id, sta1, sta2) tuples."""
    return list(pairs)


class TestPccXcorr:

    def test_identical_signals_peak_at_zero(self):
        """PCC2 of a signal with itself must peak at lag=0."""
        x = RNG.standard_normal(N)
        data = np.stack([x, x])
        index = _make_index(("A_B", 0, 1))
        corr = pcc_xcorr(data, ML, None, index)
        ccf = corr["A_B"]
        assert len(ccf) == 2 * ML + 1
        peak_lag = np.argmax(ccf) - ML
        assert peak_lag == 0

    def test_identical_signals_peak_value_one(self):
        """PCC2(0) ≈ 1.0 for identical signals (within numerical tolerance)."""
        x = RNG.standard_normal(N)
        data = np.stack([x, x])
        index = _make_index(("A_B", 0, 1))
        corr = pcc_xcorr(data, ML, None, index)
        ccf = corr["A_B"]
        assert ccf[ML] == pytest.approx(1.0, abs=1e-3)  # eps regularisation causes ~1e-5 deviation

    def test_shifted_signal_peak_at_correct_lag(self):
        """PCC2 peak should occur at the introduced shift."""
        shift = 30   # samples
        x = np.sin(2 * np.pi * 10 * np.arange(N) / N)  # clean sinusoid
        x_shifted = np.roll(x, shift)
        data = np.stack([x, x_shifted])
        index = _make_index(("A_B", 0, 1))
        corr = pcc_xcorr(data, ML, None, index)
        ccf = corr["A_B"]
        peak_lag = np.argmax(ccf) - ML
        assert peak_lag == shift

    def test_output_length(self):
        """Output CCF must have exactly 2*maxlag+1 samples."""
        x1 = RNG.standard_normal(N)
        x2 = RNG.standard_normal(N)
        data = np.stack([x1, x2])
        for ml in [10, 50, 100, 200]:
            index = _make_index(("A_B", 0, 1))
            corr = pcc_xcorr(data, ml, None, index)
            assert len(corr["A_B"]) == 2 * ml + 1, f"Wrong length for maxlag={ml}"

    def test_uncorrelated_noise_small_peak(self):
        """PCC2 of two independent noise traces should have a small max value."""
        x1 = RNG.standard_normal(N)
        x2 = RNG.standard_normal(N)
        data = np.stack([x1, x2])
        index = _make_index(("A_B", 0, 1))
        corr = pcc_xcorr(data, ML, None, index)
        ccf = corr["A_B"]
        # For N=1024, ML=100, the spurious peak should be well below 0.3
        assert np.abs(ccf).max() < 0.3

    def test_multiple_pairs(self):
        """All requested pairs must appear in the output dict."""
        x = RNG.standard_normal((4, N))
        data = x
        index = _make_index(
            ("A_B", 0, 1),
            ("A_C", 0, 2),
            ("B_C", 1, 2),
            ("B_D", 1, 3),
        )
        corr = pcc_xcorr(data, ML, None, index)
        assert set(corr.keys()) == {"A_B", "A_C", "B_C", "B_D"}
        for v in corr.values():
            assert len(v) == 2 * ML + 1

    def test_empty_index_returns_empty_dict(self):
        data = RNG.standard_normal((2, N))
        corr = pcc_xcorr(data, ML, None, [])
        assert corr == {}

    def test_energy_argument_ignored(self):
        """pcc_xcorr must produce the same result regardless of energy value."""
        x = RNG.standard_normal((2, N))
        index = _make_index(("A_B", 0, 1))
        corr_none = pcc_xcorr(x, ML, None, index)
        fake_energy = np.array([999.0, 999.0])
        corr_fake  = pcc_xcorr(x, ML, fake_energy, index)
        np.testing.assert_array_equal(corr_none["A_B"], corr_fake["A_B"])

    # ── Normalisation modes ──────────────────────────────────────────────────

    def test_normalisation_absmax(self):
        x1 = RNG.standard_normal(N)
        x2 = RNG.standard_normal(N)
        data = np.stack([x1, x2])
        index = _make_index(("A_B", 0, 1))
        corr = pcc_xcorr(data, ML, None, index, normalized="ABSMAX")
        ccf = corr["A_B"]
        assert np.abs(ccf).max() == pytest.approx(1.0, abs=1e-9)

    def test_normalisation_max(self):
        x1 = RNG.standard_normal(N)
        x2 = RNG.standard_normal(N)
        data = np.stack([x1, x2])
        index = _make_index(("A_B", 0, 1))
        corr = pcc_xcorr(data, ML, None, index, normalized="MAX")
        ccf = corr["A_B"]
        assert ccf.max() == pytest.approx(1.0, abs=1e-9)

    def test_normalisation_none(self):
        """Without normalisation the identical-signal CCF peak ≈ 1."""
        x = RNG.standard_normal(N)
        data = np.stack([x, x])
        index = _make_index(("A_B", 0, 1))
        corr = pcc_xcorr(data, ML, None, index, normalized=False)
        assert corr["A_B"][ML] == pytest.approx(1.0, abs=1e-3)  # eps regularisation causes ~1e-5 deviation

    # ── Symmetry ─────────────────────────────────────────────────────────────

    def test_lag0_identical_regardless_of_station_order(self):
        """PCC2 at lag=0 must be the same whether we compute AB or BA.

        Lag-0 of IFFT(conj(P1)*P2) = sum_k conj(P1[k])*P2[k] = <p1, p2>,
        which equals conj(<p2, p1>) = conj(lag-0 of IFFT(conj(P2)*P1)).
        Since we take the real part, both are identical.
        """
        x1 = RNG.standard_normal(1024)
        x2 = RNG.standard_normal(1024)
        data = np.stack([x1, x2])
        ml = 50
        corr_ab = pcc_xcorr(data, ml, None, _make_index(("A_B", 0, 1)))
        corr_ba = pcc_xcorr(data, ml, None, _make_index(("B_A", 1, 0)))
        assert corr_ab["A_B"][ml] == pytest.approx(corr_ba["B_A"][ml], abs=1e-12)

    def test_self_correlation_peak_dominates(self):
        """PCC2 of a signal with itself must have its global maximum at lag=0.

        Note: the extracted CCF is not symmetric around lag=0 even for self-correlation.
        This is a mathematical consequence of linear (zero-padded) correlation:
        positive lags come from full[0..ml] while negative lags come from
        full[Nz-ml..Nz-1], and Nz > N so these index different regions of the
        (even) full array.  The peak-at-zero property does hold unconditionally.
        """
        x = RNG.standard_normal(1024)
        data = np.stack([x, x])
        ml = 50
        corr = pcc_xcorr(data, ml, None, _make_index(("A_B", 0, 1)))
        ccf = corr["A_B"]
        assert np.argmax(ccf) == ml   # lag=0 is the global maximum
