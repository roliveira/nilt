"""Type stubs for the pybind11-compiled `nilt._nilt` extension module.

Every algorithm's ``__call__`` is overloaded: scalar ``t`` returns a
Python ``float``; a 1-D ``float64`` array of ``t`` returns a ``float64``
ndarray of the same length.

``Fs`` is invoked with a numpy array of s-values (real for Stehfest, complex
for Talbot / DeHoog) — including from the scalar ``t`` path, because the
pybind11 layer batches internally to keep a single round-trip. The type is
left as ``Callable[..., Any]`` because numpy broadcasting means scalar-typed
callables such as ``lambda s: 1 / (s + 1)`` work transparently.
"""

from typing import Any, Callable, overload

import numpy as np
from numpy.typing import NDArray


pi: float


class Stehfest:
    """Gaver-Stehfest algorithm (real-valued F(s))."""

    N: int
    """Number of terms (must be even, 2-20, default 18)."""

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, *, N: int = ...) -> None: ...

    @overload
    def __call__(
        self,
        Fs: Callable[..., Any],
        t: float,
    ) -> float: ...
    @overload
    def __call__(
        self,
        Fs: Callable[..., Any],
        t: NDArray[np.float64],
    ) -> NDArray[np.float64]: ...


class Talbot:
    """Fixed Talbot algorithm (complex-valued F(s))."""

    N: int
    """Number of quadrature points (default 50, table-accelerated 8-64)."""

    SHIFT: float
    """Real-axis contour shift (default 0.0)."""

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, *, N: int = ..., SHIFT: float = ...) -> None: ...

    @overload
    def __call__(
        self,
        Fs: Callable[..., Any],
        t: float,
    ) -> float: ...
    @overload
    def __call__(
        self,
        Fs: Callable[..., Any],
        t: NDArray[np.float64],
    ) -> NDArray[np.float64]: ...


class DeHoog:
    """De Hoog et al. algorithm (complex-valued F(s))."""

    M: int
    """Order of approximation (default 40)."""

    T_FACTOR: float
    """Period factor T = T_FACTOR * t (default 4.0)."""

    TOL: float
    """Contour damping tolerance (default 1e-16)."""

    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, *, M: int = ..., T_FACTOR: float = ..., TOL: float = ...) -> None: ...

    @overload
    def __call__(
        self,
        Fs: Callable[..., Any],
        t: float,
    ) -> float: ...
    @overload
    def __call__(
        self,
        Fs: Callable[..., Any],
        t: NDArray[np.float64],
    ) -> NDArray[np.float64]: ...
