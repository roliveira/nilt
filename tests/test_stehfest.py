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
        with pytest.raises(ValueError):
            nilt.invert(Fs, 0.0)

    def test_raises_for_t_negative(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, -1.0)


class TestStehfestInvalidArgument:
    
    def test_raises_for_N_odd(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, 1.0, options={"N": 5})

    def test_raises_for_N_too_small(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, 1.0, options={"N": 0})

    def test_raises_for_N_too_large(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, 1.0, options={"N": 22})


class TestStehfestDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self, Fs):
        algo = nilt.Stehfest()
        via_invert = nilt.invert(Fs, 3.0)
        via_call = algo(Fs, 3.0)
        assert via_invert == via_call


class TestStehfestIterableInput:

    def test_array_returns_ndarray_of_matching_length(self, Fs):
        result = nilt.invert(Fs, np.array([1.0, 2.0, 3.0]))
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_list_returns_ndarray(self, Fs):
        result = nilt.invert(Fs, [1.0, 2.0, 3.0])
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_tuple_returns_ndarray(self, Fs):
        result = nilt.invert(Fs, (1.0, 2.0))
        assert isinstance(result, np.ndarray)
        assert len(result) == 2

    def test_iterable_elements_match_scalar_calls(self, Fs):
        t_values = [1.0, 2.0, 3.0, 4.0]
        array_result = nilt.invert(Fs, t_values)
        for i, t in enumerate(t_values):
            scalar_result = nilt.invert(Fs, t)
            assert array_result[i] == pytest.approx(scalar_result, rel=1e-12)
