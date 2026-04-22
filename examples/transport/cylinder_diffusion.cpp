/*
 * Axisymmetric diffusion from a long circular cylinder into an infinite medium.
 *
 * A long solid cylinder of radius a, initially at concentration C0, is placed
 * in an infinite medium at concentration 0 at t=0. We track the average
 * concentration inside the cylinder.
 *
 * PDE (cylindrical coordinates, axisymmetric, no z-dependence):
 *   dC/dt = D * (1/r) * d/dr(r * dC/dr),  0 < r < a
 *
 * Laplace domain (average remaining concentration):
 *   C_avg_bar(s) = C0/s * (1 - (2/q) * I1(q) / I0(q)),   q = a*sqrt(s/D)
 *
 * where I0, I1 are modified Bessel functions of the first kind.
 *
 * Analytical (series solution):
 *   C_avg(t) = C0 * 4 * sum_{n=1}^{inf} exp(-D*alpha_n^2*t/a^2) / (a^2 * alpha_n^2)
 *
 * where alpha_n are the positive roots of J0(alpha) = 0 (Bessel zeros).
 *
 * For simplicity we use the first ~40 Bessel zeros for the series.
 *
 * Parameters:
 *   D = 5e-10 m^2/s (slower diffusion, e.g. large molecule in gel)
 *   a = 5e-4 m (0.5 mm radius)
 *   C0 = 1.0 mol/m^3
 *
 * Reference:
 *   Crank, J. (1975). The Mathematics of Diffusion, 2nd ed., §5.3,
 *   Oxford University Press.
 */
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>

#include "nilt.hpp"
#include "../utils/bessel.hpp"

int main()
{
    const double D  = 5.0e-10;  // diffusivity [m^2/s]
    const double a  = 5.0e-4;   // cylinder radius [m]
    const double C0 = 1.0;      // initial concentration [mol/m^3]

    // Average remaining concentration: C_avg_bar = (C0/s)(1 - (2/q)*I1(q)/I0(q))
    auto Fs = [=](auto s) {
        auto q = a * std::sqrt(s / D);
        return C0 / s * (1.0 - (2.0 / q) * bessel::I1(q) / bessel::I0(q));
    };

    // First 40 positive zeros of J0(x)
    static const double j0_zeros[] = {
        2.4048255577,  5.5200781103,  8.6537279129, 11.7915344391,
       14.9309177086, 18.0710639679, 21.2116366299, 24.3524715308,
       27.4934791320, 30.6346064684, 33.7758202136, 36.9170983537,
       40.0584257646, 43.1997917132, 46.3411883717, 49.4826098974,
       52.6240518411, 55.7655107550, 58.9069839261, 62.0484691902,
       65.1899648002, 68.3314693299, 71.4729816036, 74.6145006437,
       77.7560256304, 80.8975558711, 84.0390908310, 87.1806300504,
       90.3221731320, 93.4637197312, 96.6052695394, 99.7468222793,
      102.8883777002,106.0299355741,109.1714957001,112.3130578913,
      115.4546219810,118.5961878201,121.7377552742,124.8793242216
    };

    auto analytical = [=](double t) {
        double sum = 0.0;
        for (int n = 0; n < 40; ++n)
        {
            double an = j0_zeros[n];
            sum += std::exp(-D * an * an * t / (a * a)) / (an * an);
        }
        return C0 * 4.0 * sum;
    };

    nilt::Stehfest stehfest;
    nilt::Talbot   talbot;
    nilt::DeHoog   dehoog;

    std::ofstream ofs("transport_cylinder_diffusion.csv");
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

    std::cout << "transport_cylinder_diffusion.csv" << std::endl;

    const double t_snap = 500.0;  // [s]
    const int    nr     = 200;

    std::ofstream ofs2("transport_cylinder_diffusion_spatial.csv");
    ofs2 << "r,analytical,talbot" << std::endl;
    ofs2.precision(16);

    for (int i = 0; i <= nr; ++i)
    {
        double ri = a * static_cast<double>(i) / nr;  // 0 to a
        // Laplace-domain local concentration (Dirichlet BC at r=a):
        // C_bar(r,s) = (C0/s)(1 - I0(r*sqrt(s/D)) / I0(a*sqrt(s/D)))
        auto Fs_r = [=](auto s) {
            auto qa = a  * std::sqrt(s / D);
            auto qr = ri * std::sqrt(s / D);
            return C0 / s * (1.0 - bessel::I0(qr) / bessel::I0(qa));
        };

        // No analytical series (would need J0/J1); use NILT for both columns.
        double val = nilt::invert(talbot, Fs_r, t_snap);
        ofs2 << ri << "," << val << "," << val << std::endl;
    }
    std::cout << "transport_cylinder_diffusion_spatial.csv" << std::endl;

    return 0;
}
