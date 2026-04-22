/*
 * Diffusion from a sphere into an infinite medium.
 *
 * A solid sphere of radius a, initially at concentration C0, is placed in
 * an infinite medium at concentration 0 at t=0. We track the average
 * concentration inside the sphere.
 *
 * PDE (radial diffusion):
 *   dC/dt = D * (1/r^2) * d/dr(r^2 * dC/dr),  0 < r < a
 *
 * Laplace domain (average remaining concentration):
 *   C_avg_bar(s) = (C0/s) * (1 - (3/q)*(coth(q) - 1/q)),  q = a*sqrt(s/D)
 *
 * Analytical (series solution):
 *   C_avg(t) = C0 * (6/pi^2) * sum_{n=1}^{inf} (1/n^2) * exp(-D*n^2*pi^2*t/a^2)
 *
 * Parameters:
 *   D = 1e-9 m^2/s (typical ionic diffusion in water)
 *   a = 1e-3 m (1 mm radius)
 *   C0 = 1.0 mol/m^3
 *
 * Reference:
 *   Crank, J. (1975). The Mathematics of Diffusion, 2nd ed., §6.3,
 *   Oxford University Press.
 */
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>

#include "nilt.hpp"

int main()
{
    const double D  = 1.0e-9;   // diffusivity [m^2/s]
    const double a  = 1.0e-3;   // sphere radius [m]
    const double C0 = 1.0;      // initial concentration [mol/m^3]

    // Average remaining concentration:
    // C_avg_bar = (C0/s)(1 - (3/q)(coth(q) - 1/q))
    auto Fs = [=](auto s) {
        auto q = a * std::sqrt(s / D);
        auto coth_q = std::cosh(q) / std::sinh(q);
        return C0 / s * (1.0 - (3.0 / q) * (coth_q - 1.0 / q));
    };

    auto analytical = [=](double t) {
        double sum = 0.0;
        for (int n = 1; n <= 500; ++n)
        {
            double n2 = static_cast<double>(n * n);
            sum += std::exp(-D * n2 * nilt::pi * nilt::pi * t / (a * a)) / n2;
        }
        return C0 * 6.0 / (nilt::pi * nilt::pi) * sum;
    };

    nilt::Stehfest stehfest;
    nilt::Talbot   talbot;
    nilt::DeHoog   dehoog;

    std::ofstream ofs("transport_sphere_diffusion.csv");
    ofs << "t,analytical,stehfest,talbot,dehoog" << std::endl;
    ofs.precision(16);

    for (double t = 10.0; t <= 5000.0; t *= 1.06)
    {
        ofs << t
            << "," << analytical(t)
            << "," << nilt::invert(stehfest, Fs, t)
            << "," << nilt::invert(talbot,   Fs, t)
            << "," << nilt::invert(dehoog,   Fs, t)
            << std::endl;
    }

    std::cout << "transport_sphere_diffusion.csv" << std::endl;
    return 0;
}
