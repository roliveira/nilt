"""Randomized round-trip checks for the three inversion algorithms.

Uses hypothesis to draw parameters from two closed-form Laplace pairs and
verifies that ``nilt.invert`` recovers the analytical inverse to the
appropriate method-specific tolerance.

Skipped cleanly if hypothesis isn't installed so the base ``dev`` install
still passes.

The parameter ranges reflect each method's practical usability envelope:
Stehfest loses precision rapidly once ``b*t`` climbs past ~5; Talbot and
DeHoog stay accurate over a much wider range.
"""

import math

import pytest

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import assume, given, settings, strategies as st  # noqa: E402

import nilt  # noqa: E402


_HS_SETTINGS = settings(max_examples=25, deadline=None)


# pair 1: F(s) = a / (s + b),   f(t) = a * exp(-b * t)

_exp_a = st.floats(min_value=-10.0, max_value=10.0,
                   allow_nan=False, allow_infinity=False)
_exp_b = st.floats(min_value=0.1, max_value=10.0,
                   allow_nan=False, allow_infinity=False)
_exp_t = st.floats(min_value=0.1, max_value=5.0,
                   allow_nan=False, allow_infinity=False)


@given(a=_exp_a, b=_exp_b, t=_exp_t)
@_HS_SETTINGS
def test_stehfest_exponential_decay(a, b, t):
    """Stehfest loses accuracy fast once b*t > ~5; keep inside that window."""
    assume(abs(a) > 1e-3)
    assume(b * t <= 5.0)
    expected = a * math.exp(-b * t)
    got = nilt.invert(lambda s: a / (s + b), t, method="stehfest")
    assert got == pytest.approx(expected, rel=1e-3, abs=1e-6)


@given(a=_exp_a, b=_exp_b, t=_exp_t)
@_HS_SETTINGS
@pytest.mark.parametrize("algo_name,rtol,atol", [
    ("talbot", 1e-6, 1e-8),
    ("dehoog", 1e-8, 1e-10),
])
def test_talbot_dehoog_exponential_decay(algo_name, rtol, atol, a, b, t):
    assume(abs(a) > 1e-3)
    expected = a * math.exp(-b * t)
    got = nilt.invert(lambda s: a / (s + b), t, method=algo_name)
    assert got == pytest.approx(expected, rel=rtol, abs=atol)


# pair 2: F(s) = w / (s^2 + w^2),   f(t) = sin(w * t)
#
# Stehfest is skipped: oscillatory transforms need a complex-domain contour.

_sin_w = st.floats(min_value=0.5, max_value=5.0,
                   allow_nan=False, allow_infinity=False)
_sin_t = st.floats(min_value=0.1, max_value=3.0,
                   allow_nan=False, allow_infinity=False)


@given(w=_sin_w, t=_sin_t)
@_HS_SETTINGS
@pytest.mark.parametrize("algo_name,rtol,atol", [
    ("talbot", 1e-4, 1e-5),
    ("dehoog", 1e-6, 1e-8),
])
def test_talbot_dehoog_sine_pair(algo_name, rtol, atol, w, t):
    expected = math.sin(w * t)
    got = nilt.invert(lambda s: w / (s * s + w * w), t, method=algo_name)
    assert got == pytest.approx(expected, rel=rtol, abs=atol)
