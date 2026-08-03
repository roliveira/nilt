"""Top-level ``invert`` entry point."""

from nilt._registry import _select_algorithm
from nilt._validate import _prepare_t


def invert(Fs, t, method="Stehfest", **kwargs):
    """Invert the Laplace transform Fs at time(s) t.

    Parameters
    ----------
    Fs : callable
        Laplace-domain function.
        For Stehfest: ``Fs(s: float) -> float``.
        For Talbot/DeHoog: ``Fs(s: complex) -> complex``.
    t : float, int, or 1-D iterable of positive floats
        Evaluation time(s).  Must be strictly positive.  Accepts scalars,
        lists, tuples, or numpy arrays.  Complex, ``Decimal`` and
        ``Fraction`` values are rejected with ``TypeError``.
    method : str, optional
        Algorithm name: ``"Stehfest"`` (default), ``"Talbot"``, or
        ``"DeHoog"`` (case-insensitive).
    **kwargs
        Method-specific parameters.  Stehfest: ``N``.  Talbot: ``N``,
        ``SHIFT``.  DeHoog: ``M``, ``T_FACTOR``, ``TOL``.

    Returns
    -------
    float or numpy.ndarray
        ``f(t)`` — a scalar if *t* was a scalar, otherwise a 1-D
        ``float64`` array of the same length as *t*.
    """
    algo = _select_algorithm(method, kwargs)
    return algo(Fs, _prepare_t(t))
