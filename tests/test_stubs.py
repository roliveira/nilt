"""Smoke checks for the type stubs shipped with `nilt`.

These are runtime tests (not a type-checker invocation): they verify that
``py.typed`` and the ``.pyi`` files are present in the installed package
and that every symbol declared in the top-level stub actually exists at
runtime.  A separate ``tests/typing_smoke.py`` file is provided as a
mypy fixture; run it manually with ``uv run mypy tests/typing_smoke.py``.
"""

from __future__ import annotations

import ast
from pathlib import Path

import nilt


NILT_DIR = Path(nilt.__file__).parent


class TestPy_typedMarker:

    def test_py_typed_file_ships_next_to_package(self):
        assert (NILT_DIR / "py.typed").is_file()


class TestStubFilesShipped:

    def test_init_pyi_exists(self):
        assert (NILT_DIR / "__init__.pyi").is_file()

    def test_underscore_nilt_pyi_exists(self):
        assert (NILT_DIR / "_nilt.pyi").is_file()


class TestStubsMatchRuntime:
    """Every name in `__init__.pyi`'s ``__all__`` must exist at runtime."""

    def _stub_all(self) -> list[str]:
        tree = ast.parse((NILT_DIR / "__init__.pyi").read_text())
        for node in tree.body:
            if (isinstance(node, ast.Assign)
                    and len(node.targets) == 1
                    and isinstance(node.targets[0], ast.Name)
                    and node.targets[0].id == "__all__"
                    and isinstance(node.value, ast.List)):
                return [elt.value for elt in node.value.elts
                        if isinstance(elt, ast.Constant)]
        raise AssertionError("no __all__ literal in __init__.pyi")

    def test_stub_all_matches_runtime_all(self):
        assert set(self._stub_all()) == set(nilt.__all__)

    def test_every_stub_export_is_reachable_at_runtime(self):
        missing = [name for name in self._stub_all()
                   if not hasattr(nilt, name)]
        assert not missing, f"stub declares symbols missing at runtime: {missing}"
