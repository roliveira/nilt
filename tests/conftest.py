"""Shared test configuration (pytest hooks and fixtures).

Tolerance constants live in tests/tolerances.py.
This conftest ensures the tests/ directory is on sys.path so that
'from tolerances import ...' works regardless of working directory.
"""

import pytest

import os
import sys
import json
from pathlib import Path

from verification.testfunctions import TEST_FUNCTIONS


sys.path.insert(0, os.path.dirname(__file__))


@pytest.fixture(scope="session")
def TOL():
    with open(Path(__file__).parent / "tolerances.json", "r") as f:
        return json.load(f)


@pytest.fixture(scope="session")
def ft():
    return TEST_FUNCTIONS["f4"]["ft"]  # f(t) = exp(-t)


@pytest.fixture(scope="session")
def Fs():
    return TEST_FUNCTIONS["f4"]["Fs"]  # F(s) = 1/(s+1)
