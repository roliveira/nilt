import pytest
import numpy as np
import re

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

    def test_version_is_semver_ish(self):
        assert isinstance(nilt.__version__, str)
        assert nilt.__version__
        # Accept X.Y.Z with an optional pre-release / build suffix so 4.0.0rc1
        # and 0.0.0+unknown (dev fallback) both pass.
        assert re.match(r"^\d+\.\d+\.\d+", nilt.__version__)


class TestInvertRejectsInvalidMethod:

    def test_raises_for_unknown_method(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, 1.0, method="not_a_method")

    def test_lowercase_method_is_rejected(self, Fs):
        # 4.0: canonical CamelCase only (parallels scipy method strings).
        with pytest.raises(ValueError):
            nilt.invert(Fs, 1.0, method="stehfest")

    def test_unknown_option_warns_and_falls_back(self, Fs):
        # 4.0: unknown options emit UserWarning (scipy convention) instead of
        # raising - the valid keys are applied and the call still succeeds.
        with pytest.warns(UserWarning, match="Unknown option 'bad_param'"):
            result = nilt.invert(Fs, 1.0, method="Stehfest",
                                 options={"bad_param": 42})
        assert isinstance(result, float)


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

    def test_method_accepts_class(self, Fs):
        r_str = nilt.invert(Fs, 1.0, method="Talbot")
        r_cls = nilt.invert(Fs, 1.0, method=nilt.Talbot)
        assert r_str == r_cls

    def test_method_accepts_instance(self, Fs):
        r_str = nilt.invert(Fs, 1.0, method="Talbot", options={"N": 64})
        r_inst = nilt.invert(Fs, 1.0, method=nilt.Talbot(N=64))
        assert r_str == r_inst

    def test_instance_with_options_raises(self, Fs):
        with pytest.raises(TypeError, match="pre-configured"):
            nilt.invert(Fs, 1.0, method=nilt.Talbot(), options={"N": 64})

    def test_options_dict_passed_to_method(self, Fs):
        r_default = nilt.invert(Fs, 1.0)
        r_custom = nilt.invert(Fs, 1.0, options={"N": 6})
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
        r_small = nilt.invert(Fs, 1.0, options={"N": 6})
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


class TestArrayCallBatchesFs:
    """The array path must invoke Fs exactly once per call (batched)."""

    @pytest.mark.parametrize("AlgoClass", [nilt.Stehfest, nilt.Talbot, nilt.DeHoog])
    def test_array_call_invokes_Fs_once(self, AlgoClass, Fs):
        algo = AlgoClass()
        count = {"n": 0, "total": 0}

        def Fs_counted(s):
            count["n"] += 1
            count["total"] += np.asarray(s).size
            return Fs(s)

        t = np.linspace(0.5, 4.0, 25)
        algo(Fs_counted, t)

        assert count["n"] == 1
        assert count["total"] > len(t)  # batched: many s-values per t

    @pytest.mark.parametrize("AlgoClass", [nilt.Stehfest, nilt.Talbot, nilt.DeHoog])
    def test_array_call_matches_scalar_per_t(self, AlgoClass, Fs):
        algo = AlgoClass()
        t = np.array([0.25, 0.5, 1.0, 2.0, 4.0])
        arr = algo(Fs, t)
        scl = np.array([algo(Fs, float(ti)) for ti in t])
        np.testing.assert_allclose(arr, scl, rtol=1e-10, atol=1e-12)

    @pytest.mark.parametrize("AlgoClass", [nilt.Stehfest, nilt.Talbot, nilt.DeHoog])
    def test_array_call_rejects_non_positive_t(self, AlgoClass, Fs):
        algo = AlgoClass()
        t_bad = np.array([1.0, 2.0, 0.0])
        with pytest.raises((ValueError, RuntimeError)):
            algo(Fs, t_bad)


class TestMethodParamsDiscovery:
    """`_METHOD_PARAMS` mirrors the pybind11-exposed UPPER_SNAKE properties."""

    def test_stehfest_params(self):
        from nilt._registry import _METHOD_PARAMS
        assert _METHOD_PARAMS["Stehfest"] == {"N"}

    def test_talbot_params(self):
        from nilt._registry import _METHOD_PARAMS
        assert _METHOD_PARAMS["Talbot"] == {"N", "SHIFT"}

    def test_dehoog_params(self):
        from nilt._registry import _METHOD_PARAMS
        assert _METHOD_PARAMS["DeHoog"] == {"M", "T_FACTOR", "TOL"}

    def test_discovery_picks_up_new_property_on_subclass(self):
        from nilt._registry import _discover_params

        class Fake:
            EXTRA = property(lambda self: 0.0)
            N = property(lambda self: 1)
            lower_ignored = property(lambda self: None)

        assert _discover_params(Fake) == {"N", "EXTRA"}


class TestInvertBoundaryValidation:
    """A1: Python-level ``t`` validation before anything crosses into C++."""

    def test_rejects_complex_scalar(self, Fs):
        with pytest.raises(TypeError, match="complex"):
            nilt.invert(Fs, 1.0 + 0j)

    def test_rejects_numpy_complex_scalar(self, Fs):
        with pytest.raises(TypeError, match="complex"):
            nilt.invert(Fs, np.complex128(1.5))

    def test_rejects_decimal(self, Fs):
        from decimal import Decimal
        with pytest.raises(TypeError, match="Decimal"):
            nilt.invert(Fs, Decimal("1.5"))

    def test_rejects_fraction(self, Fs):
        from fractions import Fraction
        with pytest.raises(TypeError, match="Fraction"):
            nilt.invert(Fs, Fraction(3, 2))

    def test_rejects_zero_scalar(self, Fs):
        with pytest.raises(ValueError, match="positive"):
            nilt.invert(Fs, 0.0)

    def test_rejects_negative_scalar(self, Fs):
        with pytest.raises(ValueError, match="positive"):
            nilt.invert(Fs, -1.0)

    def test_rejects_2d_array(self, Fs):
        t = np.array([[1.0, 2.0], [3.0, 4.0]])
        with pytest.raises(ValueError, match="1-D"):
            nilt.invert(Fs, t)

    def test_rejects_empty_array(self, Fs):
        with pytest.raises(ValueError, match="non-empty"):
            nilt.invert(Fs, np.array([], dtype=np.float64))

    def test_rejects_array_with_zero(self, Fs):
        with pytest.raises(ValueError, match="positive"):
            nilt.invert(Fs, np.array([1.0, 0.0, 2.0]))

    def test_rejects_array_with_negative(self, Fs):
        with pytest.raises(ValueError, match="positive"):
            nilt.invert(Fs, np.array([1.0, -0.5, 2.0]))

    def test_rejects_array_with_nan(self, Fs):
        with pytest.raises(ValueError, match="finite"):
            nilt.invert(Fs, np.array([1.0, np.nan, 2.0]))

    def test_rejects_array_with_inf(self, Fs):
        with pytest.raises(ValueError, match="finite"):
            nilt.invert(Fs, np.array([1.0, np.inf, 2.0]))

    def test_accepts_python_list(self, Fs):
        got = nilt.invert(Fs, [0.5, 1.0, 2.0])
        assert isinstance(got, np.ndarray)
        assert got.shape == (3,)

    def test_accepts_python_tuple(self, Fs):
        got = nilt.invert(Fs, (0.5, 1.0, 2.0))
        assert isinstance(got, np.ndarray)
        assert got.shape == (3,)

    def test_accepts_int_scalar(self, Fs):
        got = nilt.invert(Fs, 2)
        assert isinstance(got, float)

    def test_accepts_numpy_float_scalar(self, Fs):
        got = nilt.invert(Fs, np.float64(1.5))
        assert isinstance(got, float)


class TestFsReturnValueValidation:
    """A2: ``NumpyBatched`` reports clear errors when Fs returns garbage."""

    def test_fs_returning_list_raises_type_error(self):
        # A plain Python list of floats is happily coerced by numpy but we
        # want the error surface to be predictable.  The forcecast covers
        # this case cleanly (list -> ndarray) so it succeeds; the failure
        # mode is a non-numeric return.
        algo = nilt.Stehfest()
        with pytest.raises((TypeError, ValueError)):
            algo(lambda s: "not-an-array", 1.0)

    def test_fs_returning_wrong_length_raises_value_error(self):
        algo = nilt.Stehfest()
        with pytest.raises(ValueError, match="length"):
            algo(lambda s: np.zeros(len(s) - 1), 1.0)

    def test_fs_returning_none_raises(self):
        algo = nilt.Stehfest()
        with pytest.raises((TypeError, ValueError)):
            algo(lambda s: None, 1.0)
