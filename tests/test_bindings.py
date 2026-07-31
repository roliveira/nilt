import pytest
import numpy as np

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


class TestInvertRejectsInvalidMethod:

    def test_raises_for_unknown_method(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, 1.0, method="not_a_method")

    def test_raises_for_invalid_kwarg(self, Fs):
        with pytest.raises(TypeError):
            nilt.invert(Fs, 1.0, method="Stehfest", bad_param=42)


class TestInvertScalarVsIterable:

    def test_scalar_input_returns_float(self, Fs):
        result = nilt.invert(Fs, 1.0)
        assert isinstance(result, float)

    def test_int_input_returns_float(self, Fs):
        result = nilt.invert(Fs, 1)
        assert isinstance(result, float)

    def test_ndarray_input_returns_ndarray(self, Fs):
        result = nilt.invert(Fs, np.array([1.0, 2.0, 3.0]))
        assert isinstance(result, np.ndarray)

    def test_list_input_returns_ndarray(self, Fs):
        result = nilt.invert(Fs, [1.0, 2.0, 3.0])
        assert isinstance(result, np.ndarray)

    def test_tuple_input_returns_ndarray(self, Fs):
        result = nilt.invert(Fs, (1.0, 2.0))
        assert isinstance(result, np.ndarray)

    def test_single_element_list_returns_ndarray(self, Fs):
        result = nilt.invert(Fs, [1.0])
        assert isinstance(result, np.ndarray)
        assert len(result) == 1


class TestInvertMethodSelection:

    def test_stehfest_is_default(self, Fs):
        default_result = nilt.invert(Fs, 1.0)
        explicit_result = nilt.invert(Fs, 1.0, method="Stehfest")
        assert default_result == explicit_result

    def test_method_is_case_insensitive(self, Fs):
        r1 = nilt.invert(Fs, 1.0, method="stehfest")
        r2 = nilt.invert(Fs, 1.0, method="STEHFEST")
        r3 = nilt.invert(Fs, 1.0, method="Stehfest")
        assert r1 == r2 == r3

    def test_kwargs_passed_to_method(self, Fs):
        r_default = nilt.invert(Fs, 1.0)
        r_custom = nilt.invert(Fs, 1.0, N=6)
        assert r_default != pytest.approx(r_custom, rel=1e-6)


class TestStehfestBindingExposesMutableN:

    def test_read_default_N(self):
        assert nilt.Stehfest().N == 18

    def test_write_N_persists(self):
        algo = nilt.Stehfest()
        algo.N = 10
        assert algo.N == 10

    def test_modified_N_affects_result(self, Fs):
        r_default = nilt.invert(Fs, 1.0)
        r_small = nilt.invert(Fs, 1.0, N=6)
        assert r_default != pytest.approx(r_small, rel=1e-6)


class TestTalbotBindingExposesMutableParameters:

    def test_read_default_N(self):
        assert nilt.Talbot().N == 50

    def test_read_default_SHIFT(self):
        assert nilt.Talbot().SHIFT == 0.0

    def test_write_N_persists(self):
        algo = nilt.Talbot()
        algo.N = 100
        assert algo.N == 100

    def test_write_SHIFT_persists(self):
        algo = nilt.Talbot()
        algo.SHIFT = 2.0
        assert algo.SHIFT == 2.0


class TestDeHoogBindingExposesMutableParameters:

    def test_read_default_M(self):
        assert nilt.DeHoog().M == 40

    def test_read_default_T_FACTOR(self):
        assert nilt.DeHoog().T_FACTOR == 4.0

    def test_read_default_TOL(self):
        assert nilt.DeHoog().TOL == 1e-16

    def test_write_M_persists(self):
        algo = nilt.DeHoog()
        algo.M = 20
        assert algo.M == 20

    def test_write_T_FACTOR_persists(self):
        algo = nilt.DeHoog()
        algo.T_FACTOR = 8.0
        assert algo.T_FACTOR == 8.0

    def test_write_TOL_persists(self):
        algo = nilt.DeHoog()
        algo.TOL = 1e-10
        assert algo.TOL == 1e-10


class TestAllThreeAlgorithmsCallable:
    """Each algorithm class supports direct __call__ with scalar and array."""

    @pytest.mark.parametrize("AlgoClass", [nilt.Stehfest, nilt.Talbot, nilt.DeHoog])
    def test_scalar_call_returns_float(self, AlgoClass, Fs):
        algo = AlgoClass()
        result = algo(Fs, 1.0)
        assert isinstance(result, float)

    @pytest.mark.parametrize("AlgoClass", [nilt.Stehfest, nilt.Talbot, nilt.DeHoog])
    def test_array_call_returns_ndarray(self, AlgoClass, Fs):
        algo = AlgoClass()
        t = np.array([1.0, 2.0, 3.0])
        result = algo(Fs, t)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3
