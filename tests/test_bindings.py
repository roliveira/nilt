"""Tests for the pybind11 binding layer: type exposure, dispatch, and API surface."""

import numpy as np
import pytest

import nilt


class TestModuleExportsAllSymbols:

    def test_stehfest_class_exists(self):
        assert hasattr(nilt, "Stehfest")

    def test_talbot_class_exists(self):
        assert hasattr(nilt, "Talbot")

    def test_dehoog_class_exists(self):
        assert hasattr(nilt, "DeHoog")

    def test_invert_function_exists(self):
        assert callable(nilt.invert)

    def test_pi_constant_matches_numpy(self):
        assert nilt.pi == pytest.approx(np.pi, rel=1e-15)


class TestInvertDispatchRejectsInvalidAlgorithm:

    def test_raises_for_string_algorithm(self):
        with pytest.raises((TypeError, ValueError)):
            nilt.invert("not_an_algo", lambda s: 1.0 / (s + 1.0), 1.0)

    def test_raises_for_int_algorithm(self):
        with pytest.raises((TypeError, ValueError)):
            nilt.invert(42, lambda s: 1.0 / (s + 1.0), 1.0)


class TestInvertDispatchScalarVsArray:

    def test_scalar_input_returns_float(self):
        result = nilt.invert(nilt.Talbot(), lambda s: 1.0 / (s + 1.0), 1.0)
        assert isinstance(result, float)

    def test_ndarray_input_returns_ndarray(self):
        t = np.array([1.0, 2.0])
        result = nilt.invert(nilt.Talbot(), lambda s: 1.0 / (s + 1.0), t)
        assert isinstance(result, np.ndarray)

    def test_single_element_array_returns_ndarray(self):
        t = np.array([1.0])
        result = nilt.invert(nilt.Talbot(), lambda s: 1.0 / (s + 1.0), t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 1


class TestStehfestBindingExposesMutableN:

    def test_read_default_N(self):
        assert nilt.Stehfest().N == 18

    def test_write_N_persists(self):
        algo = nilt.Stehfest()
        algo.N = 10
        assert algo.N == 10

    def test_modified_N_affects_result(self):
        Fs = lambda s: 1.0 / (s + 1.0)
        algo_default = nilt.Stehfest()
        algo_small = nilt.Stehfest()
        algo_small.N = 6
        r_default = nilt.invert(algo_default, Fs, 1.0)
        r_small = nilt.invert(algo_small, Fs, 1.0)
        assert r_default != r_small


class TestTalbotBindingExposesMutableParameters:

    def test_read_default_n(self):
        assert nilt.Talbot().n == 50

    def test_read_default_shift(self):
        assert nilt.Talbot().shift == 0.0

    def test_write_n_persists(self):
        algo = nilt.Talbot()
        algo.n = 100
        assert algo.n == 100

    def test_write_shift_persists(self):
        algo = nilt.Talbot()
        algo.shift = 2.0
        assert algo.shift == 2.0


class TestDeHoogBindingExposesMutableParameters:

    def test_read_default_M(self):
        assert nilt.DeHoog().M == 40

    def test_read_default_T_factor(self):
        assert nilt.DeHoog().T_factor == 4.0

    def test_read_default_tol(self):
        assert nilt.DeHoog().tol == 1e-16

    def test_write_M_persists(self):
        algo = nilt.DeHoog()
        algo.M = 20
        assert algo.M == 20

    def test_write_T_factor_persists(self):
        algo = nilt.DeHoog()
        algo.T_factor = 8.0
        assert algo.T_factor == 8.0

    def test_write_tol_persists(self):
        algo = nilt.DeHoog()
        algo.tol = 1e-10
        assert algo.tol == 1e-10


class TestAllThreeAlgorithmsCallable:
    """Each algorithm class supports direct __call__ with scalar and array."""

    @pytest.mark.parametrize("AlgoClass", [nilt.Stehfest, nilt.Talbot, nilt.DeHoog])
    def test_scalar_call_returns_float(self, AlgoClass):
        algo = AlgoClass()
        result = algo(lambda s: 1.0 / (s + 1.0), 1.0)
        assert isinstance(result, float)

    @pytest.mark.parametrize("AlgoClass", [nilt.Stehfest, nilt.Talbot, nilt.DeHoog])
    def test_array_call_returns_ndarray(self, AlgoClass):
        algo = AlgoClass()
        t = np.array([1.0, 2.0])
        result = algo(lambda s: 1.0 / (s + 1.0), t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 2
