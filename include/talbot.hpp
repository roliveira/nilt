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

namespace nilt {

class Talbot
{
public:
    static constexpr const char* name = "Talbot";

    int    n     = 50;     // number of quadrature points
    double shift = 0.0;    // contour shift parameter

    // Evaluate the inverse Laplace transform at time t.
    // Fs must be callable as Fs(std::complex<double>) -> std::complex<double>.
    template<typename F>
    double operator()(F&& Fs, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("Talbot: t must be positive");

        double h = 2.0 * pi / n;
        std::complex<double> ans(0.0, 0.0);

        for (int k = 0; k < n; ++k)
        {
            double theta = -pi + (k + 0.5) * h;
            double ct = std::cos(0.6407 * theta);
            double st = std::sin(0.6407 * theta);

            std::complex<double> z = shift + static_cast<double>(n) / t
                * (0.5017 * theta * ct / st + std::complex<double>(-0.6122, 0.2645 * theta));

            std::complex<double> dz = static_cast<double>(n) / t
                * (-0.5017 * 0.6407 * theta / (st * st) + 0.5017 * ct / st
                   + std::complex<double>(0.0, 0.2645));

            ans += std::exp(z * t) * Fs(z) * dz;
        }

        return (h / std::complex<double>(0.0, 2.0 * pi) * ans).real();
    }
};

} // namespace nilt

#endif // NILT_TALBOT_HEADER
