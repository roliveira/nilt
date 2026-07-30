/*
    Fixed Talbot method for numerical inversion of Laplace transforms.
    Based on the optimized contour parameters from:
    J. Abate, W. Whitt. 2006. A unified framework for numerically inverting
    Laplace transforms. INFORMS J. Comput. 18, 4, 408-421.
*/
#ifndef NILT_TALBOT_HEADER
#define NILT_TALBOT_HEADER

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "util.hpp"


namespace nilt {

namespace {

// Precomputed contour geometry for a single quadrature point.
// gamma and dgamma are the t-independent parts of z and dz:
//   z  = SHIFT + (N/t) * complex(gamma_re, gamma_im)
//   dz =         (N/t) * complex(dgamma_re, dgamma_im)
struct TalbotContourPoint {
    double gamma_re;
    double gamma_im;
    double dgamma_re;
    double dgamma_im;
};

// Flat storage: rows for N=8..64, each row has 64 slots.
// Index: (N - 8) * 64 + k
struct TalbotContourTable {
    TalbotContourPoint data[57 * 64] = {};
};

constexpr TalbotContourTable generate_talbot_table() {
    TalbotContourTable table = {};

    for (int N = 8; N <= 64; ++N) {
        double h = 2.0 * nilt::util::PI / static_cast<double>(N);
        int row_offset = (N - 8) * 64;

        for (int k = 0; k < N; ++k) {
            double theta = -nilt::util::PI + (static_cast<double>(k) + 0.5) * h;
            double arg = 0.6407 * theta;
            double ct = nilt::util::constexpr_cos(arg);
            double st = nilt::util::constexpr_sin(arg);

            // Guard: sin(arg) should never be exactly zero with the half-step
            // offset grid, but constexpr evaluation may round to zero.
            if (st == 0.0) {
                table.data[row_offset + k].gamma_re  = 0.0;
                table.data[row_offset + k].gamma_im  = 0.0;
                table.data[row_offset + k].dgamma_re = 0.0;
                table.data[row_offset + k].dgamma_im = 0.0;
                continue;
            }

            double cot = ct / st;
            double csc2 = 1.0 / (st * st);

            // gamma(theta) = 0.5017*theta*cot + complex(-0.6122, 0.2645*theta)
            table.data[row_offset + k].gamma_re  = 0.5017 * theta * cot - 0.6122;
            table.data[row_offset + k].gamma_im  = 0.2645 * theta;

            // dgamma(theta) = -0.5017*0.6407*theta*csc^2 + 0.5017*cot + complex(0, 0.2645)
            table.data[row_offset + k].dgamma_re = -0.5017 * 0.6407 * theta * csc2 + 0.5017 * cot;
            table.data[row_offset + k].dgamma_im = 0.2645;
        }
    }
    return table;
}

static constexpr TalbotContourTable TALBOT_TABLE = generate_talbot_table();

} // namespace

class Talbot
{
public:
    static constexpr const char* name = "Talbot";

    int    N     = 50;     // number of quadrature points
    double SHIFT = 0.0;    // contour shift parameter

    // Evaluate the inverse Laplace transform at time t.
    template<typename F>
    double operator()(F&& Fs, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("Talbot: t must be positive");

        if (N < 1)
            throw std::invalid_argument("Talbot: N must be >= 1");

        const double scale = static_cast<double>(N) / t;
        const double h = 2.0 * nilt::util::PI / N;

        std::complex<double> ans(0.0, 0.0);

        if (N >= 8 && N <= 64) {
            // Table path: use precomputed contour geometry
            const int row_offset = (N - 8) * 64;
            for (int k = 0; k < N; ++k)
            {
                const auto& p = TALBOT_TABLE.data[row_offset + k];
                std::complex<double> z(SHIFT + scale * p.gamma_re, scale * p.gamma_im);
                std::complex<double> dz(scale * p.dgamma_re, scale * p.dgamma_im);
                ans += std::exp(z * t) * Fs(z) * dz;
            }
        } else {
            // Runtime fallback: compute contour geometry on the fly
            for (int k = 0; k < N; ++k)
            {
                double theta = -nilt::util::PI + (k + 0.5) * h;
                double ct = std::cos(0.6407 * theta);
                double st = std::sin(0.6407 * theta);

                std::complex<double> z = SHIFT + scale
                    * (0.5017 * theta * ct / st + std::complex<double>(-0.6122, 0.2645 * theta));

                std::complex<double> dz = scale
                    * (-0.5017 * 0.6407 * theta / (st * st) + 0.5017 * ct / st
                       + std::complex<double>(0.0, 0.2645));

                ans += std::exp(z * t) * Fs(z) * dz;
            }
        }

        return (h / std::complex<double>(0.0, 2.0 * nilt::util::PI) * ans).real();
    }

    // Returns the current quadrature theta values.
    std::vector<double> get_thetas() const
    {
        if (N < 1)
            throw std::invalid_argument("Invalid N: N must be >= 1");

        double h = 2.0 * nilt::util::PI / N;
        std::vector<double> thetas(N);
        for (int k = 0; k < N; ++k) {
            thetas[k] = -nilt::util::PI + (k + 0.5) * h;
        }
        return thetas;
    }
};

} // namespace nilt

#endif // NILT_TALBOT_HEADER
