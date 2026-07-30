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

} // namespace util
} // namespace nilt

#endif // NILT_UTIL_HEADER
