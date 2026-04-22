import math

import numpy as np
import pytest

import nilt


class TestStehfestInvertExpDecay:
    """Stehfest inverts F(s)=1/(s+1) -> f(t)=exp(-t)."""

    @pytest.fixture
    def algo(self):
        return nilt.Stehfest()

    @pytest.mark.parametrize("t", [1.0, 2.0, 3.0])
    def test_relative_error_below_1e4_for_small_t(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert result == pytest.approx(math.exp(-t), rel=1e-4)

    @pytest.mark.parametrize("t", [5.0, 7.0])
    def test_relative_error_below_1e2_for_medium_t(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert result == pytest.approx(math.exp(-t), rel=1e-2)

    @pytest.mark.parametrize("t", [10.0])
    def test_relative_error_below_0_1_for_large_t(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert result == pytest.approx(math.exp(-t), rel=0.1)


class TestStehfestInvertPolynomial:
    """Stehfest inverts polynomial Laplace pairs."""

    @pytest.fixture
    def algo(self):
        return nilt.Stehfest()

    @pytest.mark.parametrize("t", [1.0, 3.0, 7.0, 10.0])
    def test_ramp_1_over_s2_equals_t(self, algo, t):
        result = nilt.invert(algo, lambda s: 1.0 / (s * s), t)
        assert result == pytest.approx(t, rel=1e-6)

    @pytest.mark.parametrize("t", [1.0, 4.0, 8.0])
    def test_cubic_1_over_s4_equals_t3_over_6(self, algo, t):
        result = nilt.invert(algo, lambda s: s ** (-4.0), t)
        assert result == pytest.approx(t ** 3 / 6.0, rel=1e-4)


class TestStehfestDefaults:

    def test_default_N_is_18(self):
        assert nilt.Stehfest().N == 18

    def test_N_is_mutable(self):
        algo = nilt.Stehfest()
        algo.N = 12
        assert algo.N == 12


class TestStehfestDomainError:

    def test_raises_for_t_zero(self):
        algo = nilt.Stehfest()
        with pytest.raises(ValueError):
            nilt.invert(algo, lambda s: 1.0 / (s + 1.0), 0.0)

    def test_raises_for_t_negative(self):
        algo = nilt.Stehfest()
        with pytest.raises(ValueError):
            nilt.invert(algo, lambda s: 1.0 / (s + 1.0), -1.0)


class TestStehfestDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self):
        algo = nilt.Stehfest()
        Fs = lambda s: 1.0 / (s + 1.0)
        via_free = nilt.invert(algo, Fs, 3.0)
        via_call = algo(Fs, 3.0)
        assert via_free == via_call


class TestStehfestArrayInput:

    def test_array_returns_ndarray_of_matching_length(self):
        algo = nilt.Stehfest()
        t = np.array([1.0, 2.0, 3.0])
        result = nilt.invert(algo, lambda s: 1.0 / (s + 1.0), t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_array_elements_match_scalar_calls(self):
        algo = nilt.Stehfest()
        Fs = lambda s: 1.0 / (s + 1.0)
        t_values = np.array([1.0, 2.0, 3.0, 4.0])
        array_result = nilt.invert(algo, Fs, t_values)
        for i, t in enumerate(t_values):
            scalar_result = nilt.invert(algo, Fs, float(t))
            assert array_result[i] == scalar_result
