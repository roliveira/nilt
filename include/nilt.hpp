/*
    nilt - Numerical Inverse Laplace Transform library
    Single-header convenience include.

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

namespace nilt {
    static constexpr double pi = 3.14159265358979323846;
} // namespace nilt

#include "stehfest.hpp"
#include "talbot.hpp"
#include "dehoog.hpp"

namespace nilt {

// Unified free function: invert(algorithm, F, t)
template<typename Algo, typename F>
double invert(const Algo& algo, F&& Fs, double t)
{
    return algo(std::forward<F>(Fs), t);
}

} // namespace nilt

#endif // NILT_HEADER
