/*
 * Theis well function - drawdown from a pumping well.
 *
 * An infinite confined aquifer, initially at head h0, is pumped at constant
 * rate Q from a fully-penetrating well at r = 0 starting at t = 0.
 * We compute the drawdown s(r,t) = h0 - h(r,t).
 *
 * Two outputs are produced:
 *   1. Drawdown vs time at a fixed observation distance (time series).
 *   2. Drawdown vs distance at a fixed snapshot time (distance profile).
 *
 * Equation (Theis, 1935 - eq. 5, US customary units):
 *   s = (114.6 Q / T) * W(u),   u = 1.87 r^2 S / (T t)
 *
 * where W(u) = E_1(u) is the exponential integral (well function).
 *
 * Laplace domain (consistent ft/day units):
 *   s_bar(r,p) = Q / (2*pi*T*p) * K0(r * sqrt(p*S/T))
 *
 * Parameters (from Theis, 1935, Fig. 1):
 *   Q = 525 gal/min    (pumping rate)
 *   T = 90,000 gal/day/ft (transmissibility)
 *   S = 0.2225          (specific yield)
 *   r = 100 ft          (observation distance, time series)
 *   t = 2 days          (snapshot time, distance profile)
 *
 * Reference:
 *   Theis, C.V. (1935). The relation between the lowering of the
 *   piezometric surface and the rate and duration of discharge of a well
 *   using ground-water storage. Trans. Am. Geophys. Union, 16(2), 519-524.
 */
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>

#include "nilt.hpp"
#include "../utils/bessel.hpp"

// Exponential integral E1(u) for real positive u.
// Uses the series for small u and the continued-fraction for large u.
double expint_e1(double u)
{
    if (u <= 0.0) return 0.0;

    if (u < 1.0)
    {
        const double euler_gamma = 0.5772156649015329;
        double sum = 0.0;
        double term = u;
        for (int n = 1; n <= 60; ++n)
        {
            sum += term / static_cast<double>(n);
            term *= -u / static_cast<double>(n + 1);
        }
        return -euler_gamma - std::log(u) + sum;
    }
    else
    {
        double a = 1.0;
        double b = u + 1.0;
        double c = 1.0e30;
        double d = 1.0 / b;
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
    // Parameters in US customary units (Theis, 1935)
    const double Q_gpm = 525.0;         // pumping rate [gal/min]
    const double T_gpd = 90000.0;       // transmissibility [gal/day/ft]
    const double S     = 0.2225;        // specific yield [-]

    // Unit conversions to consistent ft/day system
    const double gpm_to_ft3day = 1440.0 * 0.133681;  // 1 gal/min -> ft^3/day
    const double gpd_to_ft2day = 0.133681;            // 1 gal/day/ft -> ft^2/day

    const double Q = Q_gpm * gpm_to_ft3day;   // ft^3/day
    const double T = T_gpd * gpd_to_ft2day;   // ft^2/day

    // Analytical drawdown in US-customary form (eq. 5)
    auto analytical = [=](double r_ft, double t_days) {
        double u = 1.87 * r_ft * r_ft * S / (T_gpd * t_days);
        return 114.6 * Q_gpm / T_gpd * expint_e1(u);
    };

    nilt::Stehfest stehfest;
    nilt::Talbot   talbot;
    nilt::DeHoog   dehoog;

    const double r_obs = 100.0;  // observation distance [ft]

    auto Fs_time = [=](auto s) {
        return Q / (2.0 * nilt::pi * T * s)
               * bessel::K0(r_obs * std::sqrt(s * S / T));
    };

    {
        std::ofstream ofs("groundwater_theis_well_time.csv");
        ofs << "t,analytical,stehfest,talbot,dehoog" << std::endl;
        ofs.precision(16);

        for (double t = 0.01; t <= 100.0; t *= 1.06)
        {
            ofs << t
                << "," << analytical(r_obs, t)
                << "," << nilt::invert(stehfest, Fs_time, t)
                << "," << nilt::invert(talbot,   Fs_time, t)
                << "," << nilt::invert(dehoog,   Fs_time, t)
                << std::endl;
        }
        std::cout << "groundwater_theis_well_time.csv" << std::endl;
    }

    const double t_snap = 2.0;  // snapshot time [days] (48 h)

    {
        std::ofstream ofs("groundwater_theis_well_distance.csv");
        ofs << "r,analytical,stehfest,talbot,dehoog" << std::endl;
        ofs.precision(16);

        for (double r = 1.0; r <= 1000.0; r += 3.0)
        {
            auto Fs_r = [=](auto s) {
                return Q / (2.0 * nilt::pi * T * s)
                       * bessel::K0(r * std::sqrt(s * S / T));
            };
            ofs << r
                << "," << analytical(r, t_snap)
                << "," << nilt::invert(stehfest, Fs_r, t_snap)
                << "," << nilt::invert(talbot,   Fs_r, t_snap)
                << "," << nilt::invert(dehoog,   Fs_r, t_snap)
                << std::endl;
        }
        std::cout << "groundwater_theis_well_distance.csv" << std::endl;
    }

    return 0;
}
