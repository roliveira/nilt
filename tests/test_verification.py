import pytest
import numpy as np

import nilt

from verification.testfunctions import TEST_FUNCTIONS, T_VALUES


class TestVerification:
    @pytest.mark.parametrize(
        "algo_class, test_function, t",
        [
            (algo_class, key, t)
            for algo_class in [nilt.Stehfest, nilt.Talbot, nilt.DeHoog]
            for key in TEST_FUNCTIONS.keys()
            for t in T_VALUES
        ]
    )
    def test_relative_error_below_threshold(self, algo_class, test_function, t, TOL):
        algo = algo_class()
        name = algo_class.__name__
        result = nilt.invert(algo, TEST_FUNCTIONS[test_function]["Fs"], t)
        expected = TEST_FUNCTIONS[test_function]["ft"](t)

        if name == "Stehfest" and test_function in ["f4", "f7", "f8", "f9", "f10"]:
            pytest.skip("Stehfest diverges for oscillatory functions")

        assert result == pytest.approx(expected, rel=1e-4), (
            f"{name}, {test_function}(t={t}): res {result}, exp {expected}"
        )
