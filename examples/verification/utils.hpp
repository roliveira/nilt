#ifndef NILT_UTILS_HEADER
#define NILT_UTILS_HEADER

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <string>

#include "nilt.hpp"
#include "testfunctions.hpp"


template<typename Algo, typename FuncT, typename FuncS>
void output_algorithm_table(const std::string& fname, FuncT ft, FuncS Fs, const Algo& method)
{
    std::ofstream ofs(fname + "_" + method.name + ".csv");
    ofs << "t,fta,ftn,err" << std::endl;
    ofs.precision(16);

    for (double t = 1.0; t <= 10.0; t += 1.0)
    {
        double fta = ft(t);
        double ftn = nilt::invert(method, Fs, t);
        double err = std::abs(ftn - fta) / std::abs(fta);

        ofs << t << "," << fta << "," << ftn << "," << err << std::endl;
    }
}

template<typename Algo, typename FuncT, typename FuncS>
void run_function(const std::string& fname, int num, FuncT ft, FuncS Fs, const Algo& method)
{
    output_algorithm_table(fname + "_func" + std::to_string(num), ft, Fs, method);
}

void create_results_table(const std::string& fname)
{
    nilt::Stehfest stehfest;
    nilt::Talbot   talbot;
    nilt::DeHoog   dehoog;

    // Stehfest (real-valued F(s))
    run_function(fname, 1,  ft1<double>,  Fs1<double>,  stehfest);
    run_function(fname, 2,  ft2<double>,  Fs2<double>,  stehfest);
    run_function(fname, 3,  ft3<double>,  Fs3<double>,  stehfest);
    run_function(fname, 4,  ft4<double>,  Fs4<double>,  stehfest);
    run_function(fname, 5,  ft5<double>,  Fs5<double>,  stehfest);
    run_function(fname, 6,  ft6<double>,  Fs6<double>,  stehfest);
    run_function(fname, 7,  ft7<double>,  Fs7<double>,  stehfest);
    run_function(fname, 8,  ft8<double>,  Fs8<double>,  stehfest);
    run_function(fname, 9,  ft9<double>,  Fs9<double>,  stehfest);
    run_function(fname, 10, ft10<double>, Fs10<double>, stehfest);

    // Talbot and De Hoog (complex-valued F(s))
    using C = std::complex<double>;
    auto run_complex_methods = [&](int num, auto ft, auto Fs) {
        run_function(fname, num, ft, Fs, talbot);
        run_function(fname, num, ft, Fs, dehoog);
    };

    run_complex_methods(1,  ft1<double>,  Fs1<C>);
    run_complex_methods(2,  ft2<double>,  Fs2<C>);
    run_complex_methods(3,  ft3<double>,  Fs3<C>);
    run_complex_methods(4,  ft4<double>,  Fs4<C>);
    run_complex_methods(5,  ft5<double>,  Fs5<C>);
    run_complex_methods(6,  ft6<double>,  Fs6<C>);
    run_complex_methods(7,  ft7<double>,  Fs7<C>);
    run_complex_methods(8,  ft8<double>,  Fs8<C>);
    run_complex_methods(9,  ft9<double>,  Fs9<C>);
    run_complex_methods(10, ft10<double>, Fs10<C>);
}

#endif // NILT_UTILS_HEADER
