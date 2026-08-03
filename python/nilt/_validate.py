"""Validation of the ``t`` argument to :func:`nilt.invert`."""

import numbers
import numpy as np


def _prepare_t(t):
    """Validate ``t`` and return either a positive ``float`` or a 1-D
    contiguous ``float64`` array of positive finite values.

    Raises ``TypeError`` for complex, ``Decimal``, ``Fraction``, or other
    non-real numeric inputs, and ``ValueError`` for non-positive, non-finite,
    empty, or non-1-D inputs.
    """
    # Complex must be rejected before the scalar/array split because
    # ``complex`` is a numbers.Number and NumPy will happily coerce
    # complex arrays to float64 by dropping the imaginary part.
    if isinstance(t, (complex, np.complexfloating)):
        raise TypeError(f"t must be real; got complex value {t!r}")

    if isinstance(t, (int, float, np.floating, np.integer)) and not isinstance(t, bool):
        tf = float(t)
        if not (tf > 0.0):
            raise ValueError(f"t must be positive; got {tf}")
        return tf

    if isinstance(t, numbers.Number):
        # Decimal, Fraction, or any other exotic numeric type that isn't a
        # plain float/int/numpy scalar.  Coercing these silently through
        # ``np.asarray`` would lose precision or hide bugs, so reject.
        raise TypeError(
            f"t must be int, float, numpy scalar, or 1-D array-like; "
            f"got {type(t).__name__}"
        )

    try:
        t_arr = np.asarray(t, dtype=np.float64)
    except (TypeError, ValueError) as e:
        raise TypeError(
            f"t must be int, float, numpy scalar, or 1-D array-like of reals; "
            f"got {type(t).__name__}"
        ) from e

    if t_arr.ndim != 1:
        raise ValueError(
            f"t must be a 1-D array; got array with shape {t_arr.shape}"
        )
    if t_arr.size == 0:
        raise ValueError("t must be non-empty")
    if not np.all(np.isfinite(t_arr)):
        raise ValueError("t must contain only finite values")
    if np.any(t_arr <= 0.0):
        raise ValueError(
            f"all t values must be positive; smallest was {float(t_arr.min())}"
        )

    return np.ascontiguousarray(t_arr)
