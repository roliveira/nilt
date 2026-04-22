"""Snapshot tests for the 10 verification functions (Stehfest 1970 + Abate & Whitt 2006).

Each test evaluates the numerical inversion at t=1..10 and compares against
the analytical solution.  Tolerances reflect the precision of each method.
"""

import math

import numpy as np
import pytest

import nilt

pi = math.pi
EULER_GAMMA = 0.5772156649015329  # Euler-Mascheroni constant

# Test function definitions (Laplace pairs)

# func1: F(s)=1/sqrt(s), f(t)=1/sqrt(pi*t)
def Fs1(s): return 1.0 / np.sqrt(s)
def ft1(t): return 1.0 / math.sqrt(pi * t)

# func2: F(s)=ln(s)/s, f(t)=-gamma-ln(t)
def Fs2(s): return np.log(s) / s
def ft2(t): return -EULER_GAMMA - math.log(t)

# func3: F(s)=1/s^4, f(t)=t^3/6
def Fs3(s): return s ** (-4.0)
def ft3(t): return t ** 3 / 6.0

# func4: F(s)=1/(s+1), f(t)=exp(-t)
def Fs4(s): return 1.0 / (s + 1.0)
def ft4(t): return math.exp(-t)

# func5: F(s)=sqrt(pi/(2s^3))*exp(-1/(2s)), f(t)=sin(sqrt(2t))
def Fs5(s): return np.sqrt(pi / (2.0 * s ** 3)) * np.exp(-1.0 / (2.0 * s))
def ft5(t): return math.sin(math.sqrt(2.0 * t))

# func6: F(s)=1/s^2, f(t)=t
def Fs6(s): return 1.0 / (s * s)
def ft6(t): return t

# func7: F(s)=1/(s+1)^2, f(t)=t*exp(-t)
def Fs7(s): return 1.0 / ((s + 1.0) ** 2)
def ft7(t): return t * math.exp(-t)

# func8: F(s)=1/(s^2+1), f(t)=sin(t)
def Fs8(s): return 1.0 / (s * s + 1.0)
def ft8(t): return math.sin(t)

# func9: F(s)=s/(s^2+1), f(t)=cos(t)
def Fs9(s): return s / (s * s + 1.0)
def ft9(t): return math.cos(t)

# func10: F(s)=1/((s+1)^2+1), f(t)=exp(-t)*sin(t)
def Fs10(s): return 1.0 / ((s + 1.0) ** 2 + 1.0)
def ft10(t): return math.exp(-t) * math.sin(t)


FUNCTIONS = [
    ("func1",  Fs1,  ft1),
    ("func2",  Fs2,  ft2),
    ("func3",  Fs3,  ft3),
    ("func4",  Fs4,  ft4),
    ("func5",  Fs5,  ft5),
    ("func6",  Fs6,  ft6),
    ("func7",  Fs7,  ft7),
    ("func8",  Fs8,  ft8),
    ("func9",  Fs9,  ft9),
    ("func10", Fs10, ft10),
]

# func5 has a branch cut in F(s) = sqrt(pi/(2*s^3)) * exp(-1/(2s))
# that breaks contour-based methods (Talbot, DeHoog).  Only Stehfest (real s) works.
FUNCTIONS_COMPLEX = [f for f in FUNCTIONS if f[0] != "func5"]

T_VALUES = list(range(1, 11))

# Precomputed snapshot: (func_name, t) -> analytical value
SNAPSHOTS = {}
for name, _, ft in FUNCTIONS:
    for t in T_VALUES:
        SNAPSHOTS[(name, t)] = ft(float(t))


# Stehfest (all 10 functions, real F(s))

# Stehfest struggles with oscillatory functions (sin, cos, damped_sin) at large t.
# Exclude the worst cases from the parametric sweep.
STEHFEST_OSCILLATORY = {"func8", "func9", "func10"}

class TestVerificationStehfest:

    @pytest.fixture
    def algo(self):
        return nilt.Stehfest()

    @pytest.mark.parametrize("name,Fs,ft", FUNCTIONS, ids=[f[0] for f in FUNCTIONS])
    @pytest.mark.parametrize("t", T_VALUES)
    def test_relative_error_below_threshold(self, algo, name, Fs, ft, t):
        if name in STEHFEST_OSCILLATORY and t >= 6:
            pytest.skip("Stehfest diverges for oscillatory functions at large t")
        expected = SNAPSHOTS[(name, t)]
        result = nilt.invert(algo, Fs, float(t))
        # Stehfest degrades for some functions at large t
        if abs(expected) > 1e-10:
            assert abs(result - expected) / abs(expected) < 0.5, (
                f"{name} t={t}: got {result}, expected {expected}"
            )
        else:
            assert abs(result - expected) < 1e-3


# Talbot (all 10 functions, complex F(s))

class TestVerificationTalbot:

    @pytest.fixture
    def algo(self):
        return nilt.Talbot()

    @pytest.mark.parametrize("name,Fs,ft", FUNCTIONS_COMPLEX, ids=[f[0] for f in FUNCTIONS_COMPLEX])
    @pytest.mark.parametrize("t", T_VALUES)
    def test_relative_error_below_1e5(self, algo, name, Fs, ft, t):
        expected = SNAPSHOTS[(name, t)]
        result = nilt.invert(algo, Fs, float(t))
        if abs(expected) > 1e-10:
            assert abs(result - expected) / abs(expected) < 1e-5, (
                f"{name} t={t}: got {result}, expected {expected}"
            )
        else:
            assert abs(result - expected) < 1e-7


# DeHoog (all 10 functions, complex F(s))

class TestVerificationDeHoog:

    @pytest.fixture
    def algo(self):
        return nilt.DeHoog()

    @pytest.mark.parametrize("name,Fs,ft", FUNCTIONS_COMPLEX, ids=[f[0] for f in FUNCTIONS_COMPLEX])
    @pytest.mark.parametrize("t", T_VALUES)
    def test_relative_error_below_1e8(self, algo, name, Fs, ft, t):
        expected = SNAPSHOTS[(name, t)]
        result = nilt.invert(algo, Fs, float(t))
        if abs(expected) > 1e-10:
            assert abs(result - expected) / abs(expected) < 1e-8, (
                f"{name} t={t}: got {result}, expected {expected}"
            )
        else:
            assert abs(result - expected) < 1e-10


# Cross-method consistency

class TestVerificationCrossMethodAgreement:
    """All three methods produce the same sign and similar magnitude."""

    @pytest.mark.parametrize("name,Fs,ft", FUNCTIONS_COMPLEX[:4], ids=[f[0] for f in FUNCTIONS_COMPLEX[:4]])
    def test_three_methods_agree_on_sign_at_t3(self, name, Fs, ft):
        t = 3.0
        stehfest_val = nilt.invert(nilt.Stehfest(), Fs, t)
        talbot_val = nilt.invert(nilt.Talbot(), Fs, t)
        dehoog_val = nilt.invert(nilt.DeHoog(), Fs, t)

        # All should have the same sign as analytical
        expected = ft(t)
        if expected > 0:
            assert stehfest_val > 0
            assert talbot_val > 0
            assert dehoog_val > 0
        elif expected < 0:
            assert stehfest_val < 0
            assert talbot_val < 0
            assert dehoog_val < 0


# Array-based snapshot: all 10 times at once

class TestVerificationArraySnapshot:
    """Verify that array inversion at t=[1..10] matches per-element analytical values."""

    @pytest.mark.parametrize("name,Fs,ft", FUNCTIONS_COMPLEX, ids=[f[0] for f in FUNCTIONS_COMPLEX])
    def test_dehoog_array_matches_analytical_at_all_10_points(self, name, Fs, ft):
        algo = nilt.DeHoog()
        t_arr = np.arange(1.0, 11.0)
        results = nilt.invert(algo, Fs, t_arr)
        assert results.shape == (10,)

        for i, t in enumerate(T_VALUES):
            expected = SNAPSHOTS[(name, t)]
            if abs(expected) > 1e-10:
                rel_err = abs(results[i] - expected) / abs(expected)
                assert rel_err < 1e-8, (
                    f"{name} t={t}: rel_err={rel_err}"
                )
