"""Algorithm registry and tunable-parameter discovery."""

from nilt._nilt import Stehfest, Talbot, DeHoog


_METHODS = {
    "stehfest": Stehfest,
    "talbot":   Talbot,
    "dehoog":   DeHoog,
}


def _discover_params(cls):
    """Return the set of UPPER_SNAKE tunable attributes on ``cls``.

    Discovered at import time from the pybind11-exposed descriptors so a new
    knob added on the C++ side surfaces here automatically.
    """
    return {name for name in dir(cls)
            if name.isupper() and isinstance(getattr(cls, name), property)}


_METHOD_PARAMS = {key: _discover_params(cls) for key, cls in _METHODS.items()}


def _select_algorithm(method, kwargs):
    """Look up the algorithm class, instantiate it, and apply ``kwargs``.

    Raises ``ValueError`` for an unknown ``method`` and ``TypeError`` for an
    unknown kwarg.
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

    return algo
