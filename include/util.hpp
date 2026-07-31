#ifndef NILT_UTIL_HEADER
#define NILT_UTIL_HEADER

/// @file util.hpp
/// @brief Compile-time (constexpr) math utilities used by the NILT algorithms.
///
/// All functions here are constexpr-safe for C++14 and used to precompute
/// lookup tables (Stehfest coefficients, Talbot contour geometry) at compile
/// time, eliminating runtime overhead from the hot inversion loops.

#include <cstddef>


namespace nilt {
namespace util {

constexpr double PI = 3.14159265358979323846;

/// Compile-time integer power via repeated multiplication.
/// Exact for integer exponents (no log/exp rounding).
constexpr double constexpr_pow(double base, int exp) {
    double result = 1.0;
    for (int i = 0; i < exp; ++i) 
        result *= base;
    return result;
}

/// Compile-time factorial.  Valid for n in [0, 20] (20! < 2^63).
constexpr unsigned long long constexpr_factorial(int n) {
    unsigned long long x = 1;
    for (int i = 2; i <= n; ++i) 
        x *= i;
    return x;
}

/// Compile-time sine via Maclaurin series (15 terms, full double precision).
/// Input is range-reduced to [-PI, PI] before evaluation.
constexpr double constexpr_sin(double x) {
    while (x > PI)  x -= 2.0 * PI;
    while (x < -PI) x += 2.0 * PI;

    // sin(x) = x - x^3/3! + x^5/5! - ...  (Horner-style recurrence on term)
    double term = x;
    double sum = x;
    for (int n = 1; n <= 15; ++n) {
        term *= -x * x / static_cast<double>((2 * n) * (2 * n + 1));
        sum += term;
    }
    return sum;
}

/// Compile-time cosine via Maclaurin series (15 terms, full double precision).
/// Input is range-reduced to [-PI, PI] before evaluation.
constexpr double constexpr_cos(double x) {
    while (x > PI)  x -= 2.0 * PI;
    while (x < -PI) x += 2.0 * PI;

    // cos(x) = 1 - x^2/2! + x^4/4! - ...  (Horner-style recurrence on term)
    double term = 1.0;
    double sum = 1.0;
    for (int n = 1; n <= 15; ++n) {
        term *= -x * x / static_cast<double>((2 * n - 1) * (2 * n));
        sum += term;
    }
    return sum;
}

} // namespace util
} // namespace nilt

#endif // NILT_UTIL_HEADER
