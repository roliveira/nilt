/*
 * Pumping-injection well dipole - asymmetric 2D groundwater flow.
 *
 * A confined aquifer with two wells at different locations:
 *   - Pumping well at (x_p, y_p) extracting at rate Q_p
 *   - Injection well at (x_i, y_i) injecting at rate Q_i
 *
 * The net drawdown is the superposition of individual Theis solutions.
 * When Q_p != Q_i or wells are not symmetric about the origin, the
 * drawdown field has NO radial symmetry.
 *
 * Laplace domain (superposition):
 *   s_bar(x,y,p) = Q_p / (2*pi*T*p) * K0(r_p * sqrt(p/alpha))
 *                - Q_i / (2*pi*T*p) * K0(r_i * sqrt(p/alpha))
 *
 * where r_p, r_i are distances to the pumping and injection wells,
 * and the negative sign on Q_i represents head recovery from injection.
 *
 * Analytical (superposition of Theis):
 *   s(x,y,t) = Q_p/(4*pi*T) * E1(r_p^2*S/(4Tt))
 *            - Q_i/(4*pi*T) * E1(r_i^2*S/(4Tt))
 *
 * Parameters:
 *   T = 1e-3 m^2/s, S = 1e-4
 *   Well 1 (pumping):   Q_p = 0.010 m^3/s  at (-15, -5) m
 *   Well 2 (injection): Q_i = 0.007 m^3/s  at (+20, +8) m
 *
 * Reference:
 *   Bear, J. (1979). Hydraulics of Groundwater, §8.3, McGraw-Hill.
 */
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>

#include "nilt.hpp"
#include "../utils/bessel.hpp"

// Exponential integral E1(u) for real positive u.
double expint_e1(double u)
{
    if (u <= 0.0) return 0.0;
    if (u < 1.0)
    {
        const double euler_gamma = 0.5772156649015329;
        double sum = 0.0, term = u;
        for (int n = 1; n <= 60; ++n)
        {
            sum += term / static_cast<double>(n);
            term *= -u / static_cast<double>(n + 1);
        }
        return -euler_gamma - std::log(u) + sum;
    }
    else
    {
        double a = 1.0, b = u + 1.0;
        double c = 1.0e30, d = 1.0 / b;
        double result = d;
        for (int n = 1; n <= 60; ++n)
        {
            a = -static_cast<double>(n * n);
            b += 2.0;
            d = 1.0 / (b + a * d);
            c = b + a / c;
            double delta = c * d;
            result *= delta;
            if (std::abs(delta - 1.0) < 1e-15) break;
        }
        return std::exp(-u) * result;
    }
}

int main()
{
    const double T_aq  = 1.0e-3;   // transmissivity [m^2/s]
    const double S     = 1.0e-4;   // storativity
    const double alpha = T_aq / S;  // hydraulic diffusivity [m^2/s]

    // Well 1 - pumping
    const double Qp = 0.010;   // [m^3/s]
    const double xp = -15.0, yp = -5.0;
    // Well 2 - injection
    const double Qi = 0.007;   // [m^3/s]
    const double xi =  20.0, yi =  8.0;

    const double xo = 5.0, yo = 0.0;
    const double rp_obs = std::hypot(xo - xp, yo - yp);
    const double ri_obs = std::hypot(xo - xi, yo - yi);

    auto Fs_obs = [=](auto s) {
        auto sa = std::sqrt(s / alpha);
        return Qp / (2.0 * nilt::pi * T_aq * s) * bessel::K0(rp_obs * sa)
             - Qi / (2.0 * nilt::pi * T_aq * s) * bessel::K0(ri_obs * sa);
    };

    auto analytical_obs = [=](double t) {
        double up = rp_obs * rp_obs * S / (4.0 * T_aq * t);
        double ui = ri_obs * ri_obs * S / (4.0 * T_aq * t);
        return Qp / (4.0 * nilt::pi * T_aq) * expint_e1(up)
             - Qi / (4.0 * nilt::pi * T_aq) * expint_e1(ui);
    };

    nilt::Stehfest stehfest;
    nilt::Talbot   talbot;
    nilt::DeHoog   dehoog;

    {
        std::ofstream ofs("groundwater_well_dipole.csv");
        ofs << "t,analytical,stehfest,talbot,dehoog" << std::endl;
        ofs.precision(16);

        for (double t = 10.0; t <= 100000.0; t *= 1.06)
        {
            ofs << t
                << "," << analytical_obs(t)
                << "," << nilt::invert(stehfest, Fs_obs, t)
                << "," << nilt::invert(talbot,   Fs_obs, t)
                << "," << nilt::invert(dehoog,   Fs_obs, t)
                << std::endl;
        }
        std::cout << "groundwater_well_dipole.csv" << std::endl;
    }

    const double t_snap = 7200.0;  // 2 hours [s]
    const double half   = 60.0;    // half-extent [m]
    const int    nx     = 120;
    const int    ny     = 120;

    {
        std::ofstream ofs("groundwater_well_dipole_spatial.csv");
        ofs << "x,y,analytical,talbot" << std::endl;
        ofs.precision(12);

        for (int j = 0; j < ny; ++j)
        {
            double yy = -half + 2.0 * half * j / (ny - 1);
            for (int i = 0; i < nx; ++i)
            {
                double xx = -half + 2.0 * half * i / (nx - 1);
                double dp = std::hypot(xx - xp, yy - yp);
                double di = std::hypot(xx - xi, yy - yi);

                double up = dp * dp * S / (4.0 * T_aq * t_snap);
                double ui = di * di * S / (4.0 * T_aq * t_snap);
                double anal = Qp / (4.0 * nilt::pi * T_aq) * expint_e1(up)
                            - Qi / (4.0 * nilt::pi * T_aq) * expint_e1(ui);

                auto Fs_xy = [=](auto s) {
                    auto sa = std::sqrt(s / alpha);
                    return Qp / (2.0 * nilt::pi * T_aq * s) * bessel::K0(dp * sa)
                         - Qi / (2.0 * nilt::pi * T_aq * s) * bessel::K0(di * sa);
                };
                double nilt_val = nilt::invert(talbot, Fs_xy, t_snap);

                ofs << xx << "," << yy << "," << anal << "," << nilt_val << std::endl;
            }
        }
        std::cout << "groundwater_well_dipole_spatial.csv" << std::endl;
    }

    return 0;
}
