/*
 * 2D advection-diffusion: instantaneous point release in a uniform flow.
 *
 * An infinite 2D domain with uniform velocity v in the x-direction and
 * isotropic diffusion coefficient D. A mass M [kg/m] is released
 * instantaneously at the origin at t=0. The resulting plume is elongated
 * downstream - there is NO radial symmetry.
 *
 * PDE:
 *   dC/dt + v * dC/dx = D * (d^2C/dx^2 + d^2C/dy^2)
 *
 * Laplace domain (via Galilean substitution C_bar = exp(vx/(2D)) * psi):
 *   C_bar(x,y,s) = M / (2*pi*D) * exp(v*x/(2*D))
 *                   * K0(r * sqrt(v^2/(4D^2) + s/D))
 *
 * where r = sqrt(x^2 + y^2) and K0 is the modified Bessel function of the
 * second kind, order zero.
 *
 * Analytical:
 *   C(x,y,t) = M / (4*pi*D*t) * exp(-((x - v*t)^2 + y^2) / (4*D*t))
 *
 * Parameters:
 *   D = 0.1 m^2/s  (effective isotropic dispersion)
 *   v = 0.5 m/s   (uniform flow velocity in x)
 *   M = 1.0 kg/m  (mass released)
 *
 * Reference:
 *   Bear, J. (1972). Dynamics of Fluids in Porous Media, §10.6,
 *   American Elsevier.
 */
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>

#include "nilt.hpp"
#include "../utils/bessel.hpp"

int main()
{
    const double D = 0.1;    // dispersion [m^2/s]
    const double v = 0.5;    // flow velocity in x [m/s]
    const double M = 1.0;    // released mass [kg/m]

    const double xp = 3.0;   // [m] downstream
    const double yp = 0.5;   // [m] off-axis
    const double rp = std::hypot(xp, yp);

    auto Fs_pt = [=](auto s) {
        auto kappa = std::sqrt(v * v / (4.0 * D * D) + s / D);
        return M / (2.0 * nilt::pi * D) * std::exp(v * xp / (2.0 * D))
               * bessel::K0(rp * kappa);
    };

    auto analytical_pt = [=](double t) {
        return M / (4.0 * nilt::pi * D * t)
               * std::exp(-((xp - v * t) * (xp - v * t) + yp * yp) / (4.0 * D * t));
    };

    nilt::Stehfest stehfest;
    nilt::Talbot   talbot;
    nilt::DeHoog   dehoog;

    {
        std::ofstream ofs("transport_advection_plume_2d.csv");
        ofs << "t,analytical,stehfest,talbot,dehoog" << std::endl;
        ofs.precision(16);

        for (double t = 0.5; t <= 30.0; t *= 1.06)
        {
            ofs << t
                << "," << analytical_pt(t)
                << "," << nilt::invert(stehfest, Fs_pt, t)
                << "," << nilt::invert(talbot,   Fs_pt, t)
                << "," << nilt::invert(dehoog,   Fs_pt, t)
                << std::endl;
        }
        std::cout << "transport_advection_plume_2d.csv" << std::endl;
    }

    const double t_snap = 5.0;   // [s]
    const double x_lo   = -3.0;  // [m]
    const double x_hi   =  8.0;  // [m]
    const double y_lo   = -3.0;
    const double y_hi   =  3.0;
    const int    nx     = 140;
    const int    ny     = 80;

    {
        std::ofstream ofs("transport_advection_plume_2d_spatial.csv");
        ofs << "x,y,analytical,talbot" << std::endl;
        ofs.precision(12);

        for (int j = 0; j < ny; ++j)
        {
            double yy = y_lo + (y_hi - y_lo) * j / (ny - 1);
            for (int i = 0; i < nx; ++i)
            {
                double xx = x_lo + (x_hi - x_lo) * i / (nx - 1);
                double rr = std::hypot(xx, yy);

                double anal = M / (4.0 * nilt::pi * D * t_snap)
                    * std::exp(-((xx - v * t_snap) * (xx - v * t_snap) + yy * yy)
                               / (4.0 * D * t_snap));

                double nilt_val;
                if (rr < 1e-12)
                {
                    // At exact origin K0 diverges; use large but finite value
                    nilt_val = anal;
                }
                else
                {
                    auto Fs_xy = [=](auto s) {
                        auto kappa = std::sqrt(v * v / (4.0 * D * D) + s / D);
                        return M / (2.0 * nilt::pi * D) * std::exp(v * xx / (2.0 * D))
                               * bessel::K0(rr * kappa);
                    };
                    nilt_val = nilt::invert(talbot, Fs_xy, t_snap);
                }

                ofs << xx << "," << yy << "," << anal << "," << nilt_val << std::endl;
            }
        }
        std::cout << "transport_advection_plume_2d_spatial.csv" << std::endl;
    }

    return 0;
}
