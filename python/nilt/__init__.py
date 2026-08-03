"""nilt - Numerical Inverse Laplace Transform.

Three algorithms are available:

  - ``Stehfest``  - real-valued F(s), very fast.  Best default.
  - ``Talbot``    - complex-valued F(s), robust for oscillatory transforms.
  - ``DeHoog``    - complex-valued F(s), most accurate for difficult transforms.

Quick start::

    import nilt

    # scalar, default method (Stehfest)
    f_1 = nilt.invert(lambda s: 1/(s+1), 1.0)

    # iterable (list, tuple, ndarray, ...)
    f_t = nilt.invert(lambda s: 1/(s+1), [0.1, 0.5, 1.0, 2.0, 5.0])

    # different method with custom parameters (scipy-style)
    f_t = nilt.invert(lambda s: 1/(s+1), 1.0, method="Talbot", options={"N": 64})

    # or pass a pre-configured instance
    f_t = nilt.invert(lambda s: 1/(s+1), 1.0, method=nilt.Talbot(N=64))
"""

from nilt._nilt import Stehfest, Talbot, DeHoog, pi
from nilt._version import __version__
from nilt._invert import invert

__all__ = [
    "__version__",
    "pi",
    "invert",
    "Stehfest",
    "Talbot",
    "DeHoog",
]
