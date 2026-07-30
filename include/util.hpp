#ifndef NILT_UTIL_HEADER
#define NILT_UTIL_HEADER


#include <cstddef>


namespace nilt {
namespace util {

constexpr double PI = 3.14159265358979323846;

// Compile-time power function for integer exponents
constexpr double constexpr_pow(double base, int exp) {
    double result = 1.0;
    for (int i = 0; i < exp; ++i) 
        result *= base;
    return result;
}

// Compile-time factorial function
constexpr unsigned long long constexpr_factorial(int n) {
    unsigned long long x = 1;
    for (int i = 2; i <= n; ++i) 
        x *= i;
    return x;
}

// Compile-time sin via Taylor series (converges for all x, reduced to [-PI, PI])
constexpr double constexpr_sin(double x) {
    // Range reduction to [-PI, PI]
    while (x > PI)  x -= 2.0 * PI;
    while (x < -PI) x += 2.0 * PI;

    double term = x;
    double sum = x;
    for (int n = 1; n <= 15; ++n) {
        term *= -x * x / static_cast<double>((2 * n) * (2 * n + 1));
        sum += term;
    }
    return sum;
}

// Compile-time cos via Taylor series
constexpr double constexpr_cos(double x) {
    // Range reduction to [-PI, PI]
    while (x > PI)  x -= 2.0 * PI;
    while (x < -PI) x += 2.0 * PI;

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
