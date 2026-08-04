"""Algorithm registry, tunable-parameter discovery, and dispatcher."""

import warnings

from nilt._nilt import Stehfest, Talbot, DeHoog


_METHODS = {
    "Stehfest": Stehfest,
    "Talbot":   Talbot,
    "DeHoog":   DeHoog,
}


def _discover_params(cls):
    """Return the set of UPPER_SNAKE tunable attributes on ``cls``.

    Discovered at import time from the pybind11-exposed descriptors so a new
    knob added on the C++ side surfaces here automatically.
    """
    return {name for name in dir(cls)
            if name.isupper() and isinstance(getattr(cls, name), property)}


_METHOD_PARAMS = {key: _discover_params(cls) for key, cls in _METHODS.items()}


def _apply_options(algo, options, method_name):
    """Apply ``options`` to ``algo`` in place; warn on unknown keys.

    Matches scipy's convention (``scipy.optimize.minimize`` etc.):
    unknown option keys emit a ``UserWarning`` and are otherwise ignored.
    """
    valid = _METHOD_PARAMS[method_name]
    for key, value in options.items():
        if key in valid:
            setattr(algo, key, value)
        else:
            warnings.warn(
                f"Unknown option '{key}' for method {method_name}; ignored. "
                f"Valid options: {sorted(valid)}",
                UserWarning,
                stacklevel=3,
            )


def _select_algorithm(method, options):
    """Resolve ``method`` (string, class, or instance) into a configured algo.

    - String: canonical name in ``_METHODS`` (e.g. ``"Stehfest"``).  A fresh
      instance is created and ``options`` (if given) is applied.
    - Class: one of the algorithm classes.  A fresh instance is created and
      ``options`` (if given) is applied.
    - Instance: returned as-is; ``options`` must be ``None`` or empty.

    Raises ``ValueError`` for an unknown method string, ``TypeError`` for a
    method of the wrong type, or ``TypeError`` if ``options`` is combined
    with a pre-configured instance.
    """
    if isinstance(method, str):
        if method not in _METHODS:
            raise ValueError(
                f"Unknown method '{method}'. "
                f"Choose from: {', '.join(_METHODS)}"
            )
        algo = _METHODS[method]()
        if options:
            _apply_options(algo, options, method)
        return algo

    if isinstance(method, type):
        # A class object (e.g. ``method=Talbot``).
        matches = [name for name, cls in _METHODS.items() if cls is method]
        if not matches:
            raise TypeError(
                f"method class must be one of {list(_METHODS.values())}; "
                f"got {method!r}"
            )
        name = matches[0]
        algo = method()
        if options:
            _apply_options(algo, options, name)
        return algo

    # Assume instance.  Verify it's one of ours.
    matches = [name for name, cls in _METHODS.items() if isinstance(method, cls)]
    if not matches:
        raise TypeError(
            f"method must be a string, class, or instance of "
            f"{list(_METHODS.values())}; got {type(method).__name__}"
        )
    if options:
        raise TypeError(
            "options cannot be combined with a pre-configured method "
            "instance; either drop options or pass the method by name/class"
        )
    return method
