/*
 * Timing benchmark - measures wall-clock time per inversion as the main
 * tuning parameter is varied for each algorithm.
 *
 * Uses func4: F(s) = 1/(s+1), f(t) = exp(-t) as the test function (cheap to
 * evaluate so the timings reflect the algorithm cost, not the user function).
 *
 * Output: benchmark_timing.csv with columns
 *   method, param, time_us  (microseconds per single inversion at t=1)
 */
#include <chrono>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "nilt.hpp"

static constexpr int WARMUP  = 50;
static constexpr int REPEATS = 500;
static constexpr double T_EVAL = 1.0;

// Benchmark function: F(s) = 1/(s+1)
static double Fs_real(double s) { return 1.0 / (s + 1.0); }
static std::complex<double> Fs_cplx(std::complex<double> s) { return 1.0 / (s + 1.0); }

template<typename Algo, typename Fs>
double time_inversion(Algo& algo, Fs&& F)
{
    volatile double sink = 0.0;

    // warm-up
    for (int i = 0; i < WARMUP; ++i)
        sink = algo(F, T_EVAL);

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < REPEATS; ++i)
        sink = algo(F, T_EVAL);
    auto t1 = std::chrono::high_resolution_clock::now();

    (void)sink;
    double us = std::chrono::duration<double, std::micro>(t1 - t0).count() / REPEATS;
    return us;
}

int main()
{
    std::ofstream ofs("benchmark_timing.csv");
    ofs << "method,param,time_us" << std::endl;
    ofs.precision(6);

    // Stehfest: vary N (must be even)
    for (int N : {4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30})
    {
        nilt::Stehfest algo;
        algo.N = N;
        double us = time_inversion(algo, Fs_real);
        ofs << "Stehfest," << N << "," << us << std::endl;
        std::cout << "Stehfest  N=" << N  << "  " << us << " us" << std::endl;
    }

    // Talbot: vary n
    for (int n : {5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200})
    {
        nilt::Talbot algo;
        algo.n = n;
        double us = time_inversion(algo, Fs_cplx);
        ofs << "Talbot," << n << "," << us << std::endl;
        std::cout << "Talbot    n=" << n  << "  " << us << " us" << std::endl;
    }

    // DeHoog: vary M
    for (int M : {5, 10, 15, 20, 30, 40, 50, 60, 80, 100, 120, 150, 200})
    {
        nilt::DeHoog algo;
        algo.M = M;
        double us = time_inversion(algo, Fs_cplx);
        ofs << "DeHoog," << M << "," << us << std::endl;
        std::cout << "DeHoog    M=" << M  << "  " << us << " us" << std::endl;
    }

    return EXIT_SUCCESS;
}
