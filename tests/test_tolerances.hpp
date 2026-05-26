#ifndef NILT_TEST_TOLERANCES_HPP
#define NILT_TEST_TOLERANCES_HPP

#include <limits>

// Test tolerance constants for numerical Laplace inversion algorithms.
//
// These values reflect the expected numerical accuracy of each algorithm,
// not floating-point machine precision. The three methods are expected 
// to perform in the following order of accuracy:
//
//   Stehfest (lowest)  <  Talbot (middle)  <  DeHoog (highest)
//
// Stehfest uses real-valued arithmetic and a fixed number of terms, limiting
// its achievable accuracy. Talbot uses a complex contour deformation that
// converges faster. DeHoog uses a complex Padé approximation and achieves
// the highest accuracy for smooth transforms.
//
// Within each method, tolerances vary by test scenario:
// - "LARGE" variants are relaxed because all methods lose accuracy at large t
//   (the Laplace inversion becomes increasingly ill-conditioned).
// - Absolute tolerances (ABS) are used when expected values are near zero,
//   where relative error is less meaningful.
// - Per-function tolerances (e.g. EXP_SMALL vs RAMP) account for differences
//   in how well each algorithm handles specific transform characteristics
//   (oscillatory, polynomial, exponential decay, etc.).

// DeHoog method (complex Padé approximation, highest accuracy) - achieves
// ~10 digits for well-behaved transforms; ~9 digits at large t.
static constexpr double DEHOOG_REL_TOL       = 1e-10;
static constexpr double DEHOOG_REL_TOL_LARGE = 1e-9;
static constexpr double DEHOOG_ABS_TOL       = 1e-12;

// Talbot method (complex contour deformation, intermediate accuracy) - 
// achieves ~10 digits at small t but degrades to ~6 at large t due to
// sensitivity of the contour parameters.
static constexpr double TALBOT_REL_TOL       = 1e-10;
static constexpr double TALBOT_REL_TOL_LARGE = 1e-6;
static constexpr double TALBOT_ABS_TOL       = 1e-7;

// Stehfest method (real-valued, lowest accuracy) - accuracy depends
// heavily on the smoothness of the time-domain function. Polynomial/ramp
// functions invert well (~6 digits); exponentials degrade quickly with t,
// dropping to ~1-2 digits for large t or oscillatory signals.
static constexpr double STEHFEST_EXP_SMALL_REL_TOL  = 1e-4;
static constexpr double STEHFEST_EXP_MEDIUM_REL_TOL = 1e-2;
static constexpr double STEHFEST_EXP_LARGE_REL_TOL  = 0.1;
static constexpr double STEHFEST_RAMP_REL_TOL       = 1e-6;
static constexpr double STEHFEST_CUBIC_REL_TOL      = 1e-4;

// Verification suite tolerances - these are broader because the verification
// sweep covers all 10 standard test functions at t=1..10, including worst-case
// scenarios for each method.
static constexpr double TOL_STEHFEST = 0.1;
static constexpr double TOL_TALBOT   = 1e-6;
static constexpr double TOL_DEHOOG   = 1e-9;

// Array consistency - verifies that scalar and array code paths produce
// identical results. This is a roundoff-level check, so it scales with machine
// epsilon.
static constexpr double ARRAY_CONSISTENCY_TOL = 100 * std::numeric_limits<double>::epsilon();

#endif // NILT_TEST_TOLERANCES_HPP
