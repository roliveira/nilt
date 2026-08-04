/*
    De Hoog, Knight & Stokes. 1982. An improved method for numerical inversion
    of Laplace transforms. SIAM J. Sci. Stat. Comput. 3(3), 357-366.

    Implementation notes:
    - Uses accelerated quotient-difference (QD) algorithm to build a continued
      fraction from 2M+1 uniformly-spaced samples on a vertical line in the
      complex s-plane.
    - No precomputation is possible: all arithmetic depends on the function
      values Fs(gamma + i*k*PI/T) which vary with t.
    - Accuracy is controlled by M (more terms = better) and T_FACTOR (larger =
      fewer aliasing artifacts but slower convergence).
    - Works with complex-valued Laplace-domain functions.
    - Heap allocations: 5 vectors of size O(M) to O(M^2).  For large M,
      consider memory vs. accuracy tradeoff.
*/
#ifndef NILT_DEHOOG_HEADER
#define NILT_DEHOOG_HEADER

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "util.hpp"


namespace nilt {

class DeHoog
{
public:
    static constexpr const char* name = "DeHoog";

    int    M        = 40;       // order of approximation (2M+1 function evaluations)
    double T_FACTOR = 4.0;      // period factor: T = T_FACTOR * t (controls aliasing)
    double TOL      = 1.0e-16;  // Bromwich contour damping: gamma = -ln(TOL)/(2T)

    /// Apply a named option.  Returns true if `key` is known to this method,
    /// false otherwise.  Used by the string-dispatched `nilt::invert(F, t,
    /// "DeHoog", {...})` entry point to route the options map.
    bool set_option(const std::string& key, double value)
    {
        if (key == "M")        { M = static_cast<int>(value); return true; }
        if (key == "T_FACTOR") { T_FACTOR = value;             return true; }
        if (key == "TOL")      { TOL = value;                  return true; }
        return false;
    }

    /// Evaluate the inverse Laplace transform at time t.
    /// @param Fs  Laplace-domain function: Fs(complex<double>) -> complex<double>
    /// @param t   Evaluation time (must be > 0)
    /// @return    Approximation of f(t)
    template<typename F>
    double operator()(F&& Fs, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("DeHoog: t must be positive");

        int    twoM  = 2 * M;    // number of quadrature points
        int    cols  = M + 1;    // number of columns in QD table (M + 1)
        double T     = T_FACTOR * t;
        double gamma = -0.5 * std::log(TOL) / T;

        // Evaluate F(s) along the vertical line Re(s) = gamma
        // s_k = gamma + i*k*PI/T,  k = 0, 1, ..., 2M
        std::vector<std::complex<double>> Fc(twoM + 1);
        Fc[0] = 0.5 * Fs(std::complex<double>(gamma, 0.0));

        for (int i = 1; i <= twoM; ++i)
            Fc[i] = Fs(std::complex<double>(gamma, i * nilt::util::PI / T));

        // Quotient-difference (QD) algorithm - eq. (20) of De Hoog et al. 1982
        std::vector<std::complex<double>> e( (twoM + 1) * cols, 0.0);
        std::vector<std::complex<double>> q( (twoM + 1) * cols, 0.0);

        for (int i = 0; i < twoM; ++i)
            q[i * cols + 1] = Fc[i + 1] / Fc[i];

        for (int r = 1; r <= M; ++r)
        {
            for (int i = 2 * (M - r); i >= 0; --i)
                e[i * cols + r] = q[(i + 1) * cols + r] - q[i * cols + r] + e[(i + 1) * cols + (r - 1)];

            if (r < M)
            {
                for (int i = 2 * (M - r) - 1; i >= 0; --i)
                    q[i * cols + (r + 1)] = q[(i + 1) * cols + r] * e[(i + 1) * cols + r] / e[i * cols + r];
            }
        }

        // Populate d vector for continued fraction
        std::vector<std::complex<double>> d(twoM + 1);
        d[0] = Fc[0];
        for (int m = 1; m <= M; ++m)
        {
            d[2 * m - 1] = -q[0 * cols + m];
            d[2 * m]     = -e[0 * cols + m];
        }

        // Evaluate continued fraction via forward recurrence - eq. (21)
        // z = exp(i*PI*t/T) = exp(i*PI/T_FACTOR) - phase factor on the unit circle
        std::complex<double> z(std::cos(nilt::util::PI * t / T), std::sin(nilt::util::PI * t / T));

        std::vector<std::complex<double>> A(twoM + 2), B(twoM + 2);
        A[0] = 0.0;  A[1] = d[0];
        B[0] = 1.0;  B[1] = 1.0;

        for (int n = 2; n <= twoM + 1; ++n)
        {
            auto dz = d[n - 1] * z;
            A[n] = A[n - 1] + dz * A[n - 2];
            B[n] = B[n - 1] + dz * B[n - 2];
        }

        // Tail acceleration via residual correction - eqs. (23)-(24)
        // Improves convergence by estimating the tail of the continued fraction
        auto h2M = 0.5 * (1.0 + z * (d[twoM - 1] - d[twoM]));
        auto R2M = -h2M * (1.0 - std::sqrt(1.0 + z * d[twoM] / (h2M * h2M)));

        A[twoM + 1] = A[twoM] + R2M * A[twoM - 1];
        B[twoM + 1] = B[twoM] + R2M * B[twoM - 1];

        return (1.0 / T) * std::exp(gamma * t) * (A[twoM + 1] / B[twoM + 1]).real();
    }

    /// Batched evaluation: `Fs_batched(const complex* s, complex* out, int n)`
    /// fills `out[0..n-1]` with F(s[0..n-1]) in a single call.  Optimized for
    /// callers where per-callback overhead dominates (e.g. Python bindings).
    /// For pure C++ with a cheap Fs, prefer `operator()`.
    template<typename Fbatch>
    double eval_batched(Fbatch&& Fs_batched, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("DeHoog: t must be positive");

        int    twoM  = 2 * M;    // number of quadrature points
        int    cols  = M + 1;    // number of columns in QD table (M + 1)
        double T     = T_FACTOR * t;
        double gamma = -0.5 * std::log(TOL) / T;

        // Evaluate F(s) along the vertical line Re(s) = gamma
        // s_k = gamma + i*k*PI/T,  k = 0, 1, ..., 2M
        std::vector<std::complex<double>> s_vals(twoM + 1), Fc(twoM + 1);
        for (int i = 0; i <= twoM; ++i)
            s_vals[i] = std::complex<double>(gamma, i * nilt::util::PI / T);

        Fs_batched(s_vals.data(), Fc.data(), twoM + 1);
        Fc[0] *= 0.5;

        // Quotient-difference (QD) algorithm - eq. (20) of De Hoog et al. 1982
        std::vector<std::complex<double>> e( (twoM + 1) * cols, 0.0);
        std::vector<std::complex<double>> q( (twoM + 1) * cols, 0.0);

        for (int i = 0; i < twoM; ++i)
            q[i * cols + 1] = Fc[i + 1] / Fc[i];

        for (int r = 1; r <= M; ++r)
        {
            for (int i = 2 * (M - r); i >= 0; --i)
                e[i * cols + r] = q[(i + 1) * cols + r] - q[i * cols + r] + e[(i + 1) * cols + (r - 1)];

            if (r < M)
            {
                for (int i = 2 * (M - r) - 1; i >= 0; --i)
                    q[i * cols + (r + 1)] = q[(i + 1) * cols + r] * e[(i + 1) * cols + r] / e[i * cols + r];
            }
        }

        // Populate d vector for continued fraction
        std::vector<std::complex<double>> d(twoM + 1);
        d[0] = Fc[0];
        for (int m = 1; m <= M; ++m)
        {
            d[2 * m - 1] = -q[0 * cols + m];
            d[2 * m]     = -e[0 * cols + m];
        }

        // Evaluate continued fraction via forward recurrence - eq. (21)
        // z = exp(i*PI*t/T) = exp(i*PI/T_FACTOR) - phase factor on the unit circle
        std::complex<double> z(std::cos(nilt::util::PI * t / T), std::sin(nilt::util::PI * t / T));

        std::vector<std::complex<double>> A(twoM + 2), B(twoM + 2);
        A[0] = 0.0;  A[1] = d[0];
        B[0] = 1.0;  B[1] = 1.0;

        for (int n = 2; n <= twoM + 1; ++n)
        {
            auto dz = d[n - 1] * z;
            A[n] = A[n - 1] + dz * A[n - 2];
            B[n] = B[n - 1] + dz * B[n - 2];
        }

        // Tail acceleration via residual correction - eqs. (23)-(24)
        // Improves convergence by estimating the tail of the continued fraction
        auto h2M = 0.5 * (1.0 + z * (d[twoM - 1] - d[twoM]));
        auto R2M = -h2M * (1.0 - std::sqrt(1.0 + z * d[twoM] / (h2M * h2M)));

        A[twoM + 1] = A[twoM] + R2M * A[twoM - 1];
        B[twoM + 1] = B[twoM] + R2M * B[twoM - 1];

        return (1.0 / T) * std::exp(gamma * t) * (A[twoM + 1] / B[twoM + 1]).real();
    }

    /// Array-batched evaluation: invokes `Fs_batched` exactly once with all
    /// `nt * (2M+1)` s-values across every time in `t[0..nt-1]`. Ideal for
    /// bindings inverting whole arrays: one Python round-trip per call
    /// instead of nt.  Pure C++ callers with a cheap Fs should prefer
    /// `nilt::invert(Fs, t_vec)` which loops the fused scalar path.
    template<typename Fbatch>
    void eval_batched(Fbatch&& Fs_batched,
                      const double* t, double* out, int nt) const
    {
        for (int i = 0; i < nt; ++i)
            if (t[i] <= 0.0)
                throw std::domain_error("DeHoog: t must be positive");

        const int twoM = 2 * M;
        const int per  = twoM + 1;
        const size_t total = static_cast<size_t>(nt) * per;

        std::vector<std::complex<double>> s(total), fv(total);

        for (int i = 0; i < nt; ++i)
        {
            double T     = T_FACTOR * t[i];
            double gamma = -0.5 * std::log(TOL) / T;
            for (int k = 0; k <= twoM; ++k)
                s[i * per + k] = std::complex<double>(gamma, k * nilt::util::PI / T);
        }

        Fs_batched(s.data(), fv.data(), nt * per);

        for (int i = 0; i < nt; ++i)
            out[i] = combine_(t[i], fv.data() + i * per);
    }

private:
    /// Shared QD + continued-fraction combine used by both per-t
    /// `eval_batched` and the array overload.  `Fc` must point to 2M+1
    /// complex F(s_k) values with s_k = gamma + i*k*PI/T (k = 0..2M).
    /// The 0.5 factor on Fc[0] is applied internally.
    double combine_(double t, const std::complex<double>* Fc_in) const
    {
        int    twoM  = 2 * M;
        int    cols  = M + 1;
        double T     = T_FACTOR * t;
        double gamma = -0.5 * std::log(TOL) / T;

        std::vector<std::complex<double>> Fc(Fc_in, Fc_in + twoM + 1);
        Fc[0] *= 0.5;

        std::vector<std::complex<double>> e((twoM + 1) * cols, 0.0);
        std::vector<std::complex<double>> q((twoM + 1) * cols, 0.0);

        for (int i = 0; i < twoM; ++i)
            q[i * cols + 1] = Fc[i + 1] / Fc[i];

        for (int r = 1; r <= M; ++r)
        {
            for (int i = 2 * (M - r); i >= 0; --i)
                e[i * cols + r] = q[(i + 1) * cols + r] - q[i * cols + r] + e[(i + 1) * cols + (r - 1)];

            if (r < M)
            {
                for (int i = 2 * (M - r) - 1; i >= 0; --i)
                    q[i * cols + (r + 1)] = q[(i + 1) * cols + r] * e[(i + 1) * cols + r] / e[i * cols + r];
            }
        }

        std::vector<std::complex<double>> d(twoM + 1);
        d[0] = Fc[0];
        for (int m = 1; m <= M; ++m)
        {
            d[2 * m - 1] = -q[0 * cols + m];
            d[2 * m]     = -e[0 * cols + m];
        }

        std::complex<double> z(std::cos(nilt::util::PI * t / T), std::sin(nilt::util::PI * t / T));

        std::vector<std::complex<double>> A(twoM + 2), B(twoM + 2);
        A[0] = 0.0;  A[1] = d[0];
        B[0] = 1.0;  B[1] = 1.0;

        for (int n = 2; n <= twoM + 1; ++n)
        {
            auto dz = d[n - 1] * z;
            A[n] = A[n - 1] + dz * A[n - 2];
            B[n] = B[n - 1] + dz * B[n - 2];
        }

        auto h2M = 0.5 * (1.0 + z * (d[twoM - 1] - d[twoM]));
        auto R2M = -h2M * (1.0 - std::sqrt(1.0 + z * d[twoM] / (h2M * h2M)));

        A[twoM + 1] = A[twoM] + R2M * A[twoM - 1];
        B[twoM + 1] = B[twoM] + R2M * B[twoM - 1];

        return (1.0 / T) * std::exp(gamma * t) * (A[twoM + 1] / B[twoM + 1]).real();
    }
};

} // namespace nilt

#endif // NILT_DEHOOG_HEADER
