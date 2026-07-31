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
            nilt.invert(Fs, 0.0, method="Talbot")

    def test_raises_for_t_negative(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, -1.0, method="Talbot")


class TestTalbotInvalidArgument:

    def test_raises_for_N_less_than_1(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, 1.0, method="Talbot", N=0)

    def test_raises_for_N_negative(self, Fs):
        with pytest.raises(ValueError):
            nilt.invert(Fs, 1.0, method="Talbot", N=-5)


class TestTalbotTableAndFallback:

    def test_N_in_table_range_produces_finite(self, Fs):
        result = nilt.invert(Fs, 1.0, method="Talbot", N=64)
        assert np.isfinite(result)

    def test_N_outside_table_range_produces_finite(self, Fs):
        result = nilt.invert(Fs, 1.0, method="Talbot", N=100)
        assert np.isfinite(result)

    def test_table_and_fallback_agree_at_boundary(self, Fs):
        """N=64 (table) and N=65 (fallback) should give similar results."""
        r64 = nilt.invert(Fs, 1.0, method="Talbot", N=64)
        r65 = nilt.invert(Fs, 1.0, method="Talbot", N=65)
        assert r64 == pytest.approx(r65, rel=0.01)


class TestTalbotDirectCallMatchesFreeFunction:

    def test_direct_call_identical_to_invert(self, Fs):
        algo = nilt.Talbot()
        assert nilt.invert(Fs, 3.0, method="Talbot") == algo(Fs, 3.0)


class TestTalbotIterableInput:

    def test_array_returns_ndarray_of_matching_length(self, Fs):
        result = nilt.invert(Fs, [1.0, 2.0, 3.0], method="Talbot")
        assert isinstance(result, np.ndarray)
        assert len(result) == 3

    def test_iterable_elements_match_scalar_calls(self, Fs):
        t_values = [1.0, 2.0, 5.0]
        array_result = nilt.invert(Fs, t_values, method="Talbot")
        for i, t in enumerate(t_values):
            scalar_result = nilt.invert(Fs, t, method="Talbot")
            assert array_result[i] == pytest.approx(scalar_result, rel=1e-12)
