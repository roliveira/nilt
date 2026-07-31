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
            nilt.invert(Fs, 0.0, method="DeHoog")

    def test_raises_for_t_negative(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, -1.0, method="DeHoog")


class TestDeHoogDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self, Fs):
        algo = nilt.DeHoog()
        assert nilt.invert(Fs, 3.0, method="DeHoog") == algo(Fs, 3.0)


class TestDeHoogIterableInput:

    def test_array_returns_ndarray_of_matching_length(self, Fs):
        result = nilt.invert(Fs, [1.0, 2.0, 3.0], method="DeHoog")
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_iterable_elements_match_scalar_calls(self, Fs):
        t_values = [1.0, 2.0, 5.0]
        array_result = nilt.invert(Fs, t_values, method="DeHoog")
        for i, t in enumerate(t_values):
            scalar_result = nilt.invert(Fs, t, method="DeHoog")
            assert array_result[i] == pytest.approx(scalar_result, rel=1e-12)
