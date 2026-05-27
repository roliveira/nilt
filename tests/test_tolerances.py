"""Shared test tolerance constants.

These tolerances reflect the expected numerical accuracy of each inversion
algorithm, not floating-point machine precision. The three methods are
expected to perform in the following order of accuracy:

  Stehfest (lowest)  <  Talbot (middle)  <  DeHoog (highest)

Stehfest uses real-valued arithmetic and a fixed number of terms, limiting
its achievable accuracy. Talbot uses a complex contour deformation that
converges faster. DeHoog uses a complex Padé approximation and achieves
the highest accuracy for smooth transforms.

Within each method, tolerances vary by test scenario:
- "LARGE" variants are relaxed because all methods lose accuracy at large t
  (the Laplace inversion becomes increasingly ill-conditioned).
- Absolute tolerances (ABS) are used when expected values are near zero,
  where relative error is less meaningful.
- Per-function tolerances (e.g. EXP_SMALL vs RAMP) account for differences
  in how well each algorithm handles specific transform characteristics
  (oscillatory, polynomial, exponential decay, etc.).
"""

import sys

# DeHoog method (complex Padé approximation, highest accuracy) - achieves
# ~10 digits for well-behaved transforms; ~9 digits at large t.
DEHOOG_REL_TOL = 1e-10
DEHOOG_REL_TOL_LARGE = 1e-9
DEHOOG_ABS_TOL = 1e-12

# Talbot method (complex contour deformation, intermediate accuracy) -
# achieves ~10 digits at small t but degrades to ~6 at large t due to
# sensitivity of the contour parameters.
TALBOT_REL_TOL = 1e-10
TALBOT_REL_TOL_LARGE = 1e-6
TALBOT_ABS_TOL = 1e-7

# Stehfest method (real-valued, lowest accuracy) - accuracy depends
# heavily on the smoothness of the time-domain function. Polynomial/ramp
# functions invert well (~6 digits); exponentials degrade quickly with t,
# dropping to ~1-2 digits for large t or oscillatory signals.
STEHFEST_EXP_SMALL_REL_TOL = 1e-4
STEHFEST_EXP_MEDIUM_REL_TOL = 1e-2
STEHFEST_EXP_LARGE_REL_TOL = 0.1
STEHFEST_RAMP_REL_TOL = 1e-6
STEHFEST_CUBIC_REL_TOL = 1e-4

# Verification suite tolerances - these are broader because the verification
# sweep covers all 10 standard test functions at t=1..10, including
# worst-case scenarios for each method.
TOL_STEHFEST_REL = 0.5
TOL_STEHFEST_ABS = 1e-3
TOL_TALBOT_REL = 1e-5
TOL_TALBOT_ABS = 1e-7
TOL_DEHOOG_REL = 1e-8
TOL_DEHOOG_ABS = 1e-10

# Array consistency - verifies that scalar and array code paths produce
# identical results. This is a round-off level check, so it scales with
# machine epsilon.
ARRAY_CONSISTENCY_TOL = 100 * sys.float_info.epsilon
