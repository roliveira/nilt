/*
    nilt - Numerical Inverse Laplace Transform library (C++14, header-only)

    Three algorithms are provided, each suited to different use cases:

      Stehfest - Real-valued F(s) only.  Very fast (constexpr coefficients).
                 Best for smooth, monotonically decaying transforms.
      Talbot   - Complex F(s).  Moderate cost (constexpr contour table).
                 Robust for oscillatory and steep transforms.
      DeHoog   - Complex F(s).  Highest cost but most accurate for difficult
                 transforms (discontinuities, long-time behaviour).

    Usage:
        #include <nilt.hpp>

        // With a lambda
        double result = nilt::invert(nilt::Talbot{}, [](auto s){ return 1.0/(s+1.0); }, 1.0);

        // With a function pointer
        double result = nilt::invert(nilt::DeHoog{}, &my_laplace_func, 2.5);

        // With custom parameters
        nilt::Stehfest algo;
        algo.N = 12;
        double result = nilt::invert(algo, my_func, 1.0);
*/
#ifndef NILT_HEADER
#define NILT_HEADER

#include <cmath>
#include <complex>
#include <stdexcept>

#include "stehfest.hpp"
#include "talbot.hpp"
#include "dehoog.hpp"

namespace nilt {

/// Unified free function for numerical inversion of Laplace transforms.
/// @tparam Algo  Algorithm class (Stehfest, Talbot, or DeHoog)
/// @tparam F     Callable type for the Laplace-domain function
/// @param  algo  Algorithm instance (may carry tuning parameters)
/// @param  Fs    Laplace-domain function
/// @param  t     Time at which to evaluate the inverse
/// @return       Approximation of f(t)
template<typename Algo, typename F>
double invert(const Algo& algo, F&& Fs, double t)
{
    return algo(std::forward<F>(Fs), t);
}

} // namespace nilt

#endif // NILT_HEADER
