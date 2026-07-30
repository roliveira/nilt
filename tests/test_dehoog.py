import pytest
import numpy as np

import nilt


class TestDeHoogDefaults:

    def test_default_M_is_40(self):
        assert nilt.DeHoog().M == 40

    def test_default_T_FACTOR_is_4(self):
        assert nilt.DeHoog().T_FACTOR == 4.0

    def test_default_TOL_is_1e16(self):
        assert nilt.DeHoog().TOL == 1e-16

    def test_M_is_mutable(self):
        algo = nilt.DeHoog()
        algo.M = 60
        assert algo.M == 60

    def test_T_FACTOR_is_mutable(self):
        algo = nilt.DeHoog()
        algo.T_FACTOR = 2.0
        assert algo.T_FACTOR == 2.0

    def test_TOL_is_mutable(self):
        algo = nilt.DeHoog()
        algo.TOL = 1e-12
        assert algo.TOL == 1e-12


class TestDeHoogDomainError:

    def test_raises_for_t_zero(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(nilt.DeHoog(), Fs, 0.0)

    def test_raises_for_t_negative(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(nilt.DeHoog(), Fs, -1.0)


class TestDeHoogDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self, Fs):
        algo = nilt.DeHoog()
        assert nilt.invert(algo, Fs, 3.0) == algo(Fs, 3.0)


class TestDeHoogArrayInput:

    def test_array_returns_ndarray_of_matching_length(self, Fs):
        algo = nilt.DeHoog()
        t = np.array([1.0, 2.0, 3.0])
        result = nilt.invert(algo, Fs, t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_array_elements_match_scalar_calls(self, Fs):
        algo = nilt.DeHoog()
        t_values = np.array([1.0, 2.0, 5.0])
        array_result = nilt.invert(algo, Fs, t_values)
        for i, t in enumerate(t_values):
            scalar_result = nilt.invert(algo, Fs, float(t))
            assert array_result[i] == pytest.approx(scalar_result, rel=1e-12)
