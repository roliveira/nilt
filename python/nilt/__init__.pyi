"""Type stubs for the top-level `nilt` package."""

from typing import Any, Callable, Sequence, overload

import numpy as np
from numpy.typing import ArrayLike, NDArray

from nilt._nilt import DeHoog, Stehfest, Talbot, pi


__version__: str


@overload
def invert(
    Fs: Callable[..., Any],
    t: float | int,
    method: str = ...,
    **kwargs: Any,
) -> float: ...
@overload
def invert(
    Fs: Callable[..., Any],
    t: Sequence[float] | NDArray[np.float64] | ArrayLike,
    method: str = ...,
    **kwargs: Any,
) -> NDArray[np.float64]: ...


__all__ = [
    "Stehfest",
    "Talbot",
    "DeHoog",
    "invert",
    "pi",
    "__version__",
]
