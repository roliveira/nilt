"""Package version discovery."""

try:
    from importlib.metadata import PackageNotFoundError, version as _pkg_version
    __version__ = _pkg_version("nilt-python")
    del _pkg_version, PackageNotFoundError
except PackageNotFoundError:  # pragma: no cover - only hit in odd dev setups
    __version__ = "0.0.0+unknown"
