/*
    Harald Stehfest. 1970. Algorithm 368: Numerical inversion of Laplace transforms [D5]. 
    Commun. ACM 13, 1 (Jan. 1970), 47-49. DOI:https://doi.org/10.1145/361953.361969
*/
#ifndef NILT_STEHFEST_HEADER
#define NILT_STEHFEST_HEADER

#include <cmath>
#include <stdexcept>
#include <array>
#include <vector>
#include <algorithm>

#include "util.hpp"

namespace nilt {

namespace {

struct RawRow    { double data[21]  = {0.0}; };
struct RawMatrix { double data[210] = {0.0}; };

constexpr RawRow generate_constexpr_row(size_t N) {
    RawRow row;
    if (N < 2 || N > 20 || N % 2 != 0) 
        return row;

    size_t N2 = N / 2;
    int sign = (N2 % 2 != 0) ? -1 : 1;

    for (size_t i = 0; i < N; ++i) {
        size_t kmin = (i + 2) / 2;
        size_t kmax = std::min(i + 1, N2);
        double sum = 0.0;
        sign = -sign;

        for (size_t k = kmin; k <= kmax; ++k) {
            sum += nilt::util::constexpr_pow(static_cast<double>(k), N2)
                 * nilt::util::constexpr_factorial(2 * k)
                 / (nilt::util::constexpr_factorial(k) * nilt::util::constexpr_factorial(k - 1)
                    * nilt::util::constexpr_factorial(N2 - k) * nilt::util::constexpr_factorial(i + 1 - k)
                    * nilt::util::constexpr_factorial(2 * k - i - 1));
        }
        row.data[i] = sign * sum;
    }
    return row;
}

constexpr RawMatrix generate_constexpr_table() {
    RawMatrix matrix;
    for (size_t n = 2; n <= 20; n += 2) {
        RawRow row = generate_constexpr_row(n);
        size_t offset = (n / 2 - 1) * 21;
        for (size_t i = 0; i < 21; ++i) {
            matrix.data[offset + i] = row.data[i];
        }
    }
    return matrix;
}

static constexpr RawMatrix CONST_COEFFICIENT_RAW_MATRIX = generate_constexpr_table();

} // namespace

class Stehfest
{
public:
    static constexpr const char* name = "Stehfest";
    
    int N = 18;  // number of terms (must be even)

    // Evaluate the inverse Laplace transform at time t.
    template<typename F>
    double operator()(F&& Fs, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("Stehfest: t must be positive");
        
        if (N < 2 || N > 20 || N % 2 != 0) 
            throw std::invalid_argument("Invalid N: N must be even and between 2 and 20");

        // Reference the array directly from the helper matrix
        const double* coeff = &CONST_COEFFICIENT_RAW_MATRIX.data[(N / 2 - 1) * 21];
        double ln2t = std::log(2.0) / t;
        double s = 0.0;
        double y = 0.0;

        for (int i = 0; i < N; ++i)
        {
            s += ln2t;
            y += coeff[i] * Fs(s);
        }

        return ln2t * y;
    }

    // Returns a vector containing the N coefficients currently in use.
    std::vector<double> get_coefficients() const
    {
        if (N < 2 || N > 20 || N % 2 != 0) 
            throw std::invalid_argument("Invalid N: N must be even and between 2 and 20");

        const double* coeff_ptr = &CONST_COEFFICIENT_RAW_MATRIX.data[(N / 2 - 1) * 21];
        
        // Copy only the N valid elements into a vector
        return std::vector<double>(coeff_ptr, coeff_ptr + N);
    }
};

} // namespace nilt

#endif // NILT_STEHFEST_HEADER
