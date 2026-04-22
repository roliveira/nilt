import math

import numpy as np
import pytest

import nilt


class TestDeHoogInvertExpDecay:
    """DeHoog inverts F(s)=1/(s+1) -> f(t)=exp(-t)."""

    @pytest.fixture
    def algo(self):
        return nilt.DeHoog()

    @pytest.mark.parametrize("t", [1.0, 2.0, 5.0])
    def test_relative_error_below_1e12(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert result == pytest.approx(math.exp(-t), rel=1e-12)

    @pytest.mark.parametrize("t", [10.0])
    def test_relative_error_below_1e9_for_large_t(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert result == pytest.approx(math.exp(-t), rel=1e-9)


class TestDeHoogInvertSinusoids:

    @pytest.fixture
    def algo(self):
        return nilt.DeHoog()

    @pytest.mark.parametrize("t", [1.0, 2.0, 5.0, 9.0])
    def test_sin_t_within_1e12(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s * s + 1.0), t)
        assert result == pytest.approx(math.sin(t), abs=1e-12)

    @pytest.mark.parametrize("t", [1.0, 3.0, 6.0])
    def test_cos_t_within_1e12(self, algo, t):
        result = nilt.invert(algo, lambda s: s / (s * s + 1.0), t)
        assert result == pytest.approx(math.cos(t), abs=1e-12)

    @pytest.mark.parametrize("t", [1.0, 4.0, 8.0])
    def test_damped_sin_exp_neg_t_sin_t_within_1e12(self, algo, t):
        result = nilt.invert(
            algo, lambda s: 1.0 / ((s + 1.0) ** 2 + 1.0), t
        )
        expected = math.exp(-t) * math.sin(t)
        assert result == pytest.approx(expected, abs=1e-12)


class TestDeHoogInvertPolynomial:

    @pytest.fixture
    def algo(self):
        return nilt.DeHoog()

    @pytest.mark.parametrize("t", [1.0, 3.0, 7.0, 10.0])
    def test_ramp_1_over_s2_equals_t_within_1e12(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s * s), t)
        assert result == pytest.approx(t, rel=1e-12)


class TestDeHoogDefaults:

    def test_default_M_is_40(self):
        assert nilt.DeHoog().M == 40

    def test_default_T_factor_is_4(self):
        assert nilt.DeHoog().T_factor == 4.0

    def test_default_tol_is_1e16(self):
        assert nilt.DeHoog().tol == 1e-16

    def test_M_is_mutable(self):
        algo = nilt.DeHoog()
        algo.M = 60
        assert algo.M == 60

    def test_T_factor_is_mutable(self):
        algo = nilt.DeHoog()
        algo.T_factor = 2.0
        assert algo.T_factor == 2.0

    def test_tol_is_mutable(self):
        algo = nilt.DeHoog()
        algo.tol = 1e-12
        assert algo.tol == 1e-12


class TestDeHoogDomainError:

    def test_raises_for_t_zero(self):
        with pytest.raises(ValueError):
            nilt.invert(nilt.DeHoog(), lambda s: 1.0 / (s + 1.0), 0.0)

    def test_raises_for_t_negative(self):
        with pytest.raises(ValueError):
            nilt.invert(nilt.DeHoog(), lambda s: 1.0 / (s + 1.0), -1.0)


class TestDeHoogDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self):
        algo = nilt.DeHoog()
        Fs = lambda s: 1.0 / (s + 1.0)
        assert nilt.invert(algo, Fs, 3.0) == algo(Fs, 3.0)


class TestDeHoogArrayInput:

    def test_array_returns_ndarray_of_matching_length(self):
        algo = nilt.DeHoog()
        t = np.array([1.0, 2.0, 3.0])
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_array_elements_match_scalar_calls(self):
        algo = nilt.DeHoog()
        Fs = lambda s: 1.0 / (s + 1.0)
        t_values = np.array([1.0, 2.0, 5.0])
        array_result = nilt.invert(algo, Fs, t_values)
        for i, t in enumerate(t_values):
            assert array_result[i] == nilt.invert(algo, Fs, float(t))
