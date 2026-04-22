import math

import numpy as np
import pytest

import nilt


class TestTalbotInvertExpDecay:
    """Talbot inverts F(s)=1/(s+1) -> f(t)=exp(-t)."""

    @pytest.fixture
    def algo(self):
        return nilt.Talbot()

    @pytest.mark.parametrize("t", [1.0, 2.0, 3.0])
    def test_relative_error_below_1e10_for_small_t(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert result == pytest.approx(math.exp(-t), rel=1e-10)

    @pytest.mark.parametrize("t", [5.0, 10.0])
    def test_relative_error_below_1e6_for_large_t(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert result == pytest.approx(math.exp(-t), rel=1e-6)


class TestTalbotInvertSinusoids:

    @pytest.fixture
    def algo(self):
        return nilt.Talbot()

    @pytest.mark.parametrize("t", [1.0, 2.0, 5.0, 9.0])
    def test_sin_t_within_1e7(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s * s + 1.0), t)
        assert result == pytest.approx(math.sin(t), abs=1e-7)

    @pytest.mark.parametrize("t", [1.0, 3.0, 6.0])
    def test_cos_t_within_1e7(self, algo, t):
        result = nilt.invert(algo, lambda s: s / (s * s + 1.0), t)
        assert result == pytest.approx(math.cos(t), abs=1e-7)

    @pytest.mark.parametrize("t", [1.0, 4.0, 8.0])
    def test_damped_sin_exp_neg_t_sin_t_within_1e7(self, algo, t):
        result = nilt.invert(
            algo, lambda s: 1.0 / ((s + 1.0) ** 2 + 1.0), t
        )
        expected = math.exp(-t) * math.sin(t)
        assert result == pytest.approx(expected, abs=1e-7)


class TestTalbotInvertPolynomial:

    @pytest.fixture
    def algo(self):
        return nilt.Talbot()

    @pytest.mark.parametrize("t", [1.0, 3.0, 7.0, 10.0])
    def test_ramp_1_over_s2_equals_t_within_1e10(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s * s), t)
        assert result == pytest.approx(t, rel=1e-10)


class TestTalbotDefaults:

    def test_default_n_is_50(self):
        assert nilt.Talbot().n == 50

    def test_default_shift_is_zero(self):
        assert nilt.Talbot().shift == 0.0

    def test_n_is_mutable(self):
        algo = nilt.Talbot()
        algo.n = 100
        assert algo.n == 100

    def test_shift_is_mutable(self):
        algo = nilt.Talbot()
        algo.shift = 1.5
        assert algo.shift == 1.5


class TestTalbotDomainError:

    def test_raises_for_t_zero(self):
        with pytest.raises(ValueError):
            nilt.invert(nilt.Talbot(), lambda s: 1.0 / (s + 1.0), 0.0)

    def test_raises_for_t_negative(self):
        with pytest.raises(ValueError):
            nilt.invert(nilt.Talbot(), lambda s: 1.0 / (s + 1.0), -1.0)


class TestTalbotDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self):
        algo = nilt.Talbot()
        Fs = lambda s: 1.0 / (s + 1.0)
        assert nilt.invert(algo, Fs, 3.0) == algo(Fs, 3.0)


class TestTalbotArrayInput:

    def test_array_returns_ndarray_of_matching_length(self):
        algo = nilt.Talbot()
        t = np.array([1.0, 2.0, 3.0])
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_array_elements_match_scalar_calls(self):
        algo = nilt.Talbot()
        Fs = lambda s: 1.0 / (s + 1.0)
        t_values = np.array([1.0, 2.0, 5.0])
        array_result = nilt.invert(algo, Fs, t_values)
        for i, t in enumerate(t_values):
            assert array_result[i] == nilt.invert(algo, Fs, float(t))
