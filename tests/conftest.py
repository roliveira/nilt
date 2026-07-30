"""Shared pytest fixtures for nilt unit tests."""

import pytest

from verification.testfunctions import TEST_FUNCTIONS


@pytest.fixture(scope="session")
def ft():
    return TEST_FUNCTIONS["f4"]["ft"]  # f(t) = exp(-t)


@pytest.fixture(scope="session")
def Fs():
    return TEST_FUNCTIONS["f4"]["Fs"]  # F(s) = 1/(s+1)
