import pytest
import numpy as np

import nilt


class TestStehfestDefaults:

    def test_default_N_is_18(self):
        assert nilt.Stehfest().N == 18

    def test_N_is_mutable(self):
        algo = nilt.Stehfest()
        algo.N = 12
        assert algo.N == 12


class TestStehfestDomainError:

    def test_raises_for_t_zero(self, Fs):
        algo = nilt.Stehfest()
        with pytest.raises(ValueError):
            nilt.invert(algo, Fs, 0.0)

    def test_raises_for_t_negative(self, Fs):
        algo = nilt.Stehfest()
        with pytest.raises(ValueError):
            nilt.invert(algo, Fs, -1.0)


class TestStehfestInvalidArgument:
    
    def test_raises_for_N_odd(self, Fs):
        algo = nilt.Stehfest()
        algo.N = 5
        with pytest.raises(ValueError):
            nilt.invert(algo, Fs, 1.0)

    def test_raises_for_N_too_small(self, Fs):
        algo = nilt.Stehfest()
        algo.N = 0
        with pytest.raises(ValueError):
            nilt.invert(algo, Fs, 1.0)

    def test_raises_for_N_too_large(self, Fs):
        algo = nilt.Stehfest()
        algo.N = 22
        with pytest.raises(ValueError):
            nilt.invert(algo, Fs, 1.0)


class TestStehfestDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self, Fs):
        algo = nilt.Stehfest()
        via_free = nilt.invert(algo, Fs, 3.0)
        via_call = algo(Fs, 3.0)
        assert via_free == via_call


class TestStehfestArrayInput:

    def test_array_returns_ndarray_of_matching_length(self, Fs):
        algo = nilt.Stehfest()
        t = np.array([1.0, 2.0, 3.0])
        result = nilt.invert(algo, Fs, t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_array_elements_match_scalar_calls(self, Fs):
        algo = nilt.Stehfest()
        t_values = np.array([1.0, 2.0, 3.0, 4.0])
        array_result = nilt.invert(algo, Fs, t_values)
        for i, t in enumerate(t_values):
            scalar_result = nilt.invert(algo, Fs, float(t))
            assert array_result[i] == pytest.approx(scalar_result, rel=1e-12)
