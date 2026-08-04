"""Type stubs for the top-level `nilt` package."""

from typing import Any, Callable, Mapping, Sequence, Type, Union, overload

import numpy as np
from numpy.typing import ArrayLike, NDArray

from nilt._nilt import DeHoog, Stehfest, Talbot, pi


__version__: str

_AlgoT = Union[Stehfest, Talbot, DeHoog]
_MethodT = Union[str, Type[_AlgoT], _AlgoT]


@overload
def invert(
    Fs: Callable[..., Any],
    t: float | int,
    method: _MethodT = ...,
    options: Mapping[str, float] | None = ...,
) -> float: ...
@overload
def invert(
    Fs: Callable[..., Any],
    t: Sequence[float] | NDArray[np.float64] | ArrayLike,
    method: _MethodT = ...,
    options: Mapping[str, float] | None = ...,
) -> NDArray[np.float64]: ...


__all__ = [
    "Stehfest",
    "Talbot",
    "DeHoog",
    "invert",
    "pi",
    "__version__",
]
