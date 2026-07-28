import pytest
import numpy as np

import nilt


class TestTalbotDefaults:

    def test_default_N_is_50(self):
        assert nilt.Talbot().N == 50

    def test_default_SHIFT_is_zero(self):
        assert nilt.Talbot().SHIFT == 0.0

    def test_N_is_mutable(self):
        algo = nilt.Talbot()
        algo.N = 100
        assert algo.N == 100

    def test_SHIFT_is_mutable(self):
        algo = nilt.Talbot()
        algo.SHIFT = 1.5
        assert algo.SHIFT == 1.5


class TestTalbotDomainError:

    def test_raises_for_t_zero(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(nilt.Talbot(), Fs, 0.0)

    def test_raises_for_t_negative(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(nilt.Talbot(), Fs, -1.0)


class TestTalbotDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self, Fs):
        algo = nilt.Talbot()
        assert nilt.invert(algo, Fs, 3.0) == algo(Fs, 3.0)


class TestTalbotArrayInput:

    def test_array_returns_ndarray_of_matching_length(self, Fs):
        algo = nilt.Talbot()
        t = np.array([1.0, 2.0, 3.0])
        result = nilt.invert(algo, Fs, t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_array_elements_match_scalar_calls(self, TOL, Fs):
        algo = nilt.Talbot()
        t_values = np.array([1.0, 2.0, 5.0])
        array_result = nilt.invert(algo, Fs, t_values)
        for i, t in enumerate(t_values):
            scalar_result = nilt.invert(algo, Fs, float(t))
            assert array_result[i] == pytest.approx(scalar_result, rel=TOL["TALBOT_REL_TOL_LARGE"])
