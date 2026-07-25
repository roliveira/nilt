"""Snapshot tests for numerical inverse Laplace transform methods.

These tests lock in the exact numerical output (at full double precision) for
all 10 standard verification functions x 3 methods.  They complement the
tolerance-based tests: tolerance tests verify analytical correctness, snapshots
detect any unintended change in numerical output.

Update snapshots after intentional algorithm changes:
    pytest --snapshot-update tests/test_snapshots.py
"""

import math
from pathlib import Path

import numpy as np
import pytest

import nilt

pi = math.pi
EULER_GAMMA = 0.5772156649015329

# Laplace-domain functions
def Fs1(s): return 1.0 / np.sqrt(s)
def Fs2(s): return np.log(s) / s
def Fs3(s): return s ** (-4.0)
def Fs4(s): return 1.0 / (s + 1.0)
def Fs5(s): return np.sqrt(pi / (2.0 * s**3)) * np.exp(-1.0 / (2.0 * s))
def Fs6(s): return 1.0 / (s * s)
def Fs7(s): return 1.0 / ((s + 1.0) ** 2)
def Fs8(s): return 1.0 / (s * s + 1.0)
def Fs9(s): return s / (s * s + 1.0)
def Fs10(s): return 1.0 / ((s + 1.0) ** 2 + 1.0)

FUNCTIONS = [
    ("func1", Fs1),
    ("func2", Fs2),
    ("func3", Fs3),
    ("func4", Fs4),
    ("func5", Fs5),
    ("func6", Fs6),
    ("func7", Fs7),
    ("func8", Fs8),
    ("func9", Fs9),
    ("func10", Fs10),
]

# func5 has a branch cut that breaks contour-based methods
FUNCTIONS_COMPLEX = [(n, f) for n, f in FUNCTIONS if n != "func5"]

T_VALUES = [float(t) for t in range(1, 11)]

SNAPSHOT_DIR = Path(__file__).parent / "snapshots"


def _serialize(method, Fs):
    """Compute inversions at t=1..10 and serialize to CSV string."""
    lines = ["t,ftn"]
    for t in T_VALUES:
        ftn = nilt.invert(method, Fs, t)
        lines.append(f"{t:.1f},{ftn:.8e}")
    return "\n".join(lines) + "\n"


class TestSnapshotStehfest:

    @pytest.fixture
    def algo(self):
        return nilt.Stehfest()

    @pytest.mark.parametrize("name,Fs", FUNCTIONS, ids=[f[0] for f in FUNCTIONS])
    def test_snapshot(self, algo, name, Fs, snapshot):
        snapshot.snapshot_dir = str(SNAPSHOT_DIR)
        result = _serialize(algo, Fs)
        snapshot.assert_match(result, f"Stehfest_{name}.txt")


class TestSnapshotTalbot:

    @pytest.fixture
    def algo(self):
        return nilt.Talbot()

    @pytest.mark.parametrize("name,Fs", FUNCTIONS_COMPLEX, ids=[f[0] for f in FUNCTIONS_COMPLEX])
    def test_snapshot(self, algo, name, Fs, snapshot):
        snapshot.snapshot_dir = str(SNAPSHOT_DIR)
        result = _serialize(algo, Fs)
        snapshot.assert_match(result, f"Talbot_{name}.txt")


class TestSnapshotDeHoog:

    @pytest.fixture
    def algo(self):
        return nilt.DeHoog()

    @pytest.mark.parametrize("name,Fs", FUNCTIONS_COMPLEX, ids=[f[0] for f in FUNCTIONS_COMPLEX])
    def test_snapshot(self, algo, name, Fs, snapshot):
        snapshot.snapshot_dir = str(SNAPSHOT_DIR)
        result = _serialize(algo, Fs)
        snapshot.assert_match(result, f"DeHoog_{name}.txt")
