"""nilt - Numerical Inverse Laplace Transform.

Three algorithms are available:

  - ``Stehfest``  - real-valued F(s), very fast.  Best default.
  - ``Talbot``    - complex-valued F(s), robust for oscillatory transforms.
  - ``DeHoog``    - complex-valued F(s), most accurate for difficult transforms.

Quick start::

    import nilt

    # scalar
    f_1 = nilt.invert(lambda s: 1/(s+1), 1.0)

    # iterable (list, tuple, ndarray, ...)
    f_t = nilt.invert(lambda s: 1/(s+1), [0.1, 0.5, 1.0, 2.0, 5.0])

    # different method with custom parameters
    f_t = nilt.invert(lambda s: 1/(s+1), 1.0, method="Talbot", N=64)
"""

import numbers
import numpy as np

from nilt._nilt import Stehfest, Talbot, DeHoog, pi

try:
    from importlib.metadata import PackageNotFoundError, version as _pkg_version
    __version__ = _pkg_version("nilt-python")
    del _pkg_version, PackageNotFoundError
except PackageNotFoundError:  # pragma: no cover - only hit in odd dev setups
    __version__ = "0.0.0+unknown"


_METHODS = {
    "stehfest": Stehfest,
    "talbot":   Talbot,
    "dehoog":   DeHoog,
}

_METHOD_PARAMS = {
    "stehfest": {"N"},
    "talbot":   {"N", "SHIFT"},
    "dehoog":   {"M", "T_FACTOR", "TOL"},
}


def invert(Fs, t, method="Stehfest", **kwargs):
    """Invert the Laplace transform Fs at time(s) t.

    Parameters
    ----------
    Fs : callable
        Laplace-domain function.
        For Stehfest: ``Fs(s: float) -> float``.
        For Talbot/DeHoog: ``Fs(s: complex) -> complex``.
    t : float, int, or iterable of floats
        Evaluation time(s).  Must be positive.  Accepts scalars, lists,
        tuples, numpy arrays, or any iterable that ``numpy.asarray`` can
        convert to a 1-D float64 array.
    method : str, optional
        Algorithm name: ``"Stehfest"`` (default), ``"Talbot"``, or
        ``"DeHoog"`` (case-insensitive).
    **kwargs
        Method-specific parameters.  Stehfest: ``N``.  Talbot: ``N``,
        ``SHIFT``.  DeHoog: ``M``, ``T_FACTOR``, ``TOL``.

    Returns
    -------
    float or numpy.ndarray
        ``f(t)`` — a scalar if *t* was a scalar, otherwise an array.
    """
    key = method.lower()
    if key not in _METHODS:
        raise ValueError(
            f"Unknown method '{method}'. "
            f"Choose from: {', '.join(m.title() for m in _METHODS)}"
        )

    algo = _METHODS[key]()

    valid = _METHOD_PARAMS[key]
    for k, v in kwargs.items():
        if k not in valid:
            raise TypeError(
                f"'{k}' is not a valid parameter for {method}. "
                f"Valid parameters: {', '.join(sorted(valid))}"
            )
        setattr(algo, k, v)

    if isinstance(t, numbers.Number):
        return algo(Fs, float(t))

    t_arr = np.ascontiguousarray(t, dtype=np.float64).ravel()
    return algo(Fs, t_arr)


__all__ = [
    "Stehfest",
    "Talbot",
    "DeHoog",
    "invert",
    "pi",
    "__version__",
]
