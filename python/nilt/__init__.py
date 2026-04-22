"""nilt - Numerical Inverse Laplace Transform.

Provides three algorithms for numerically inverting Laplace transforms:

  - ``Stehfest``  - real-valued F(s), Gaver–Stehfest method.
  - ``Talbot``    - complex-valued F(s), fixed Talbot contour.
  - ``DeHoog``    - complex-valued F(s), De Hoog et al. accelerated Fourier.

Quick start::

    import numpy as np
    import nilt

    # scalar
    f_1 = nilt.invert(nilt.Talbot(), lambda s: 1/(s+1), 1.0)

    # array
    t = np.linspace(0.1, 5.0, 100)
    f_t = nilt.invert(nilt.Talbot(), lambda s: 1/(s+1), t)
"""

from nilt._nilt import Stehfest, Talbot, DeHoog, invert, pi


__all__ = [
    "Stehfest", 
    "Talbot", 
    "DeHoog", 
    "invert", 
    "pi",
]
