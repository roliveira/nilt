/*
    Harald Stehfest. 1970. Algorithm 368: Numerical inversion of Laplace transforms [D5]. 
    Commun. ACM 13, 1 (Jan. 1970), 47-49. DOI:https://doi.org/10.1145/361953.361969

    Implementation notes:
    - Coefficients for all valid N (2,4,...,20) are precomputed at compile time
      in a flat constexpr array (CONST_COEFFICIENT_RAW_MATRIX).  This avoids the
      O(N^2) factorial/power computation at runtime - the hot loop reduces to a
      single pointer dereference per coefficient (~120x faster at N=18).
    - The algorithm only works for real-valued Laplace-domain functions.
    - Accuracy degrades for large t due to exponential spacing of abscissae;
      N=18 is a good default for most problems.
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

// Storage for one row of Stehfest coefficients (max N=20, padded to 21).
struct RawRow    { double data[21]  = {0.0}; };

// Flat storage for all 10 coefficient rows: N=2,4,...,20.
// Layout: row index = (N/2 - 1), each row occupies 21 doubles.
// Total: 10 rows * 21 = 210 doubles.
struct RawMatrix { double data[210] = {0.0}; };

/// Compute Stehfest coefficients for a given even N at compile time.
/// Formula: V_i = (-1)^(N/2+i) * sum_{k=floor((i+1)/2)}^{min(i,N/2)}
///              k^(N/2) * (2k)! / (k! * (k-1)! * (N/2-k)! * (i-k)! * (2k-i-1)!)
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
    
    int N = 18;  // number of terms (must be even, 2 <= N <= 20)

    /// Evaluate the inverse Laplace transform at time t.
    /// @param Fs  Laplace-domain function: Fs(double s) -> double
    /// @param t   Evaluation time (must be > 0)
    /// @return    Approximation of f(t)
    template<typename F>
    double operator()(F&& Fs, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("Stehfest: t must be positive");

        if (N < 2 || N > 20 || N % 2 != 0)
            throw std::invalid_argument("Invalid N: N must be even and between 2 and 20");

        // Direct pointer into the precomputed coefficient table - no allocation
        const double* coeff = &CONST_COEFFICIENT_RAW_MATRIX.data[(N / 2 - 1) * 21];
        double ln2t = std::log(2.0) / t;
        double s = 0.0;
        double y = 0.0;

        // f(t) ~ (ln2/t) * sum_{i=1}^{N} V_i * F(i * ln2/t)
        for (int i = 0; i < N; ++i)
        {
            s += ln2t;
            y += coeff[i] * Fs(s);
        }

        return ln2t * y;
    }

    /// Batched evaluation: `Fs_batched(const double* s, double* out, int n)`
    /// fills `out[0..n-1]` with F(s[0..n-1]) in a single call.  Optimized for
    /// callers where per-callback overhead dominates (e.g. Python bindings):
    /// reduces N scalar callbacks to one vectorized call.  For pure C++ with
    /// a cheap inlinable Fs, prefer `operator()` - it is significantly faster
    /// because the fused loop keeps intermediates in registers.
    template<typename Fbatch>
    double eval_batched(Fbatch&& Fs_batched, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("Stehfest: t must be positive");

        if (N < 2 || N > 20 || N % 2 != 0)
            throw std::invalid_argument("Invalid N: N must be even and between 2 and 20");

        const double* coeff = &CONST_COEFFICIENT_RAW_MATRIX.data[(N / 2 - 1) * 21];
        double ln2t = std::log(2.0) / t;

        double s_vals[20];
        double f_vals[20];
        for (int i = 0; i < N; ++i)
            s_vals[i] = (i + 1) * ln2t;

        Fs_batched(s_vals, f_vals, N);

        double y = 0.0;
        for (int i = 0; i < N; ++i)
            y += coeff[i] * f_vals[i];

        return ln2t * y;
    }

    /// Array-batched evaluation: invokes `Fs_batched` exactly once with all
    /// `nt * N` s-values across every time in `t[0..nt-1]`.  Ideal for
    /// bindings inverting whole arrays: one Python round-trip per call
    /// instead of nt.  Pure C++ callers with a cheap Fs should prefer
    /// `nilt::invert(algo, Fs, t_vec)` which loops the fused scalar path.
    template<typename Fbatch>
    void eval_batched(Fbatch&& Fs_batched,
                      const double* t, double* out, int nt) const
    {
        if (N < 2 || N > 20 || N % 2 != 0)
            throw std::invalid_argument("Invalid N: N must be even and between 2 and 20");
        for (int i = 0; i < nt; ++i)
            if (t[i] <= 0.0)
                throw std::domain_error("Stehfest: t must be positive");

        std::vector<double> s(static_cast<size_t>(nt) * N);
        std::vector<double> fv(static_cast<size_t>(nt) * N);

        for (int i = 0; i < nt; ++i)
        {
            double ln2t = std::log(2.0) / t[i];
            for (int k = 0; k < N; ++k)
                s[i * N + k] = (k + 1) * ln2t;
        }

        Fs_batched(s.data(), fv.data(), nt * N);

        const double* coeff = &CONST_COEFFICIENT_RAW_MATRIX.data[(N / 2 - 1) * 21];
        for (int i = 0; i < nt; ++i)
        {
            double ln2t = std::log(2.0) / t[i];
            double y = 0.0;
            for (int k = 0; k < N; ++k)
                y += coeff[k] * fv[i * N + k];
            out[i] = ln2t * y;
        }
    }

    /// Returns a copy of the N precomputed Stehfest coefficients.
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
