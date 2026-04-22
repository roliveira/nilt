/*
 * bessel.hpp - Modified Bessel functions for use in nilt examples.
 *
 * Provides I0(z), I1(z), K0(z) for complex arguments, needed by the
 * 2D/axisymmetric examples (Theis well, cylinder diffusion).
 *
 * I0, I1: power series (converges for all z, fast for moderate |z|).
 * K0: connection formula K0(z) = -(ln(z/2)+gamma)*I0(z) + S(z) for |z|<=8,
 *      asymptotic expansion for |z|>8.
 */
#ifndef NILT_EXAMPLES_BESSEL_HPP
#define NILT_EXAMPLES_BESSEL_HPP

#include <cmath>
#include <complex>

namespace bessel {

static constexpr double euler_gamma = 0.5772156649015329;

// Modified Bessel function I0(z) via power series.
template<typename T>
T I0(T z)
{
    T sum  = 1.0;
    T term = 1.0;
    T z2_4 = z * z * 0.25;
    for (int k = 1; k <= 50; ++k)
    {
        term *= z2_4 / (static_cast<double>(k) * static_cast<double>(k));
        sum += term;
        if (std::abs(term) < std::abs(sum) * 1e-16) break;
    }
    return sum;
}

// Modified Bessel function I1(z) via power series.
template<typename T>
T I1(T z)
{
    T sum  = 1.0;
    T term = 1.0;
    T z2_4 = z * z * 0.25;
    for (int k = 1; k <= 50; ++k)
    {
        term *= z2_4 / (static_cast<double>(k) * static_cast<double>(k + 1));
        sum += term;
        if (std::abs(term) < std::abs(sum) * 1e-16) break;
    }
    return sum * z * 0.5;
}

// Modified Bessel function K0(z) for complex z.
// Uses the series representation for |z| <= 8, asymptotic for |z| > 8.
template<typename T>
T K0(T z)
{
    if (std::abs(z) <= 8.0)
    {
        // K0(z) = -(ln(z/2) + gamma) * I0(z) + sum_{k=1}^inf (z/2)^{2k} * H_k / (k!)^2
        // where H_k = 1 + 1/2 + ... + 1/k (harmonic numbers).
        T ln_z2 = std::log(z * 0.5);
        T i0 = I0(z);

        T z2_4 = z * z * 0.25;
        T term = 1.0;   // (z/2)^{2k} / (k!)^2, starts at k=0 = 1
        double Hk = 0.0;
        T series = 0.0; // sum starts at k=1
        for (int k = 1; k <= 50; ++k)
        {
            term *= z2_4 / (static_cast<double>(k) * static_cast<double>(k));
            Hk += 1.0 / static_cast<double>(k);
            series += term * Hk;
            if (std::abs(term * Hk) < std::abs(series) * 1e-16) break;
        }
        return -(ln_z2 + euler_gamma) * i0 + series;
    }
    else
    {
        // Asymptotic expansion: K0(z) ~ sqrt(pi/(2z)) * exp(-z) * sum
        // Coefficients: a_k = product_{m=1}^{k} (4*0-(2m-1)^2) / (8m)
        // i.e. a_k = product_{m=1}^k (-(2m-1)^2/(8m))
        T iz = 1.0 / z;
        T sum = 1.0;
        T ak  = 1.0;
        for (int k = 1; k <= 12; ++k)
        {
            double odd = 2.0 * k - 1.0;
            ak *= -(odd * odd) / (8.0 * static_cast<double>(k));
            T contrib = ak;
            for (int j = 0; j < k; ++j) contrib *= iz;
            sum += contrib;
        }
        return std::sqrt(nilt::pi / (2.0 * z)) * std::exp(-z) * sum;
    }
}

} // namespace bessel

#endif // NILT_EXAMPLES_BESSEL_HPP
