"""Shared test configuration (pytest hooks and fixtures).

Tolerance constants live in tests/tolerances.py.
This conftest ensures the tests/ directory is on sys.path so that
'from tolerances import ...' works regardless of working directory.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
