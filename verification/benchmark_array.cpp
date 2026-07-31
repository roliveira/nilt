/*
 * Array-size timing benchmark - measures wall-clock time for whole-array
 * inversions as len(t) is varied.  Algorithm parameters are kept at their
 * defaults (Stehfest N=18, Talbot N=50, DeHoog M=40).
 *
 * Test function: F(s) = 1/(s+1),  f(t) = exp(-t).
 *
 * Output: cpp_benchmark_array.csv with columns
 *   method, nt, time_us   (microseconds per whole-array inversion)
 */
#include <chrono>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

#include "nilt.hpp"

static constexpr int WARMUP  = 20;
static constexpr int REPEATS = 100;

static double Fs_real(double s) { return 1.0 / (s + 1.0); }
static std::complex<double> Fs_cplx(std::complex<double> s) { return 1.0 / (s + 1.0); }

template<typename Algo, typename Fs>
double time_array(const Algo& algo, Fs&& F, const std::vector<double>& t)
{
    volatile double sink = 0.0;

    for (int i = 0; i < WARMUP; ++i) {
        auto r = nilt::invert(algo, F, t);
        sink = r.empty() ? 0.0 : r.front();
    }

    auto t0 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < REPEATS; ++i) {
        auto r = nilt::invert(algo, F, t);
        sink = r.empty() ? 0.0 : r.front();
    }
    auto t1 = std::chrono::high_resolution_clock::now();

    (void)sink;
    return std::chrono::duration<double, std::micro>(t1 - t0).count() / REPEATS;
}

static std::vector<double> make_t(int nt)
{
    std::vector<double> t(nt);
    for (int i = 0; i < nt; ++i)
        t[i] = 0.1 + (5.0 - 0.1) * (nt == 1 ? 0.5 : double(i) / (nt - 1));
    return t;
}

int main()
{
    const std::vector<int> sizes = {1, 5, 10, 50, 100, 500, 1000, 5000, 10000};

    std::ofstream ofs("cpp_benchmark_array.csv");
    ofs << "method,nt,time_us" << std::endl;
    ofs.precision(6);

    nilt::Stehfest steh;
    nilt::Talbot   tal;
    nilt::DeHoog   dh;

    for (int nt : sizes) {
        auto t = make_t(nt);
        double us = time_array(steh, Fs_real, t);
        ofs << "Stehfest," << nt << "," << us << std::endl;
        std::cout << "Stehfest  nt=" << nt << "  " << us << " us" << std::endl;
    }
    for (int nt : sizes) {
        auto t = make_t(nt);
        double us = time_array(tal, Fs_cplx, t);
        ofs << "Talbot," << nt << "," << us << std::endl;
        std::cout << "Talbot    nt=" << nt << "  " << us << " us" << std::endl;
    }
    for (int nt : sizes) {
        auto t = make_t(nt);
        double us = time_array(dh, Fs_cplx, t);
        ofs << "DeHoog," << nt << "," << us << std::endl;
        std::cout << "DeHoog    nt=" << nt << "  " << us << " us" << std::endl;
    }

    return EXIT_SUCCESS;
}
