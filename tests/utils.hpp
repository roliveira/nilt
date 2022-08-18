#ifndef NILT_UTILS_HEADER
#define NILT_UTILS_HEADER

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "invlap.hpp"
#include "stehfest.hpp"
#include "dehoog.hpp"
#include "talbot.hpp"

#include "testfunctions.hpp"


void output_algorithm_table(std::string fname, double (*ft)(double), double (*Fs)(double), InverseLaplaceAlgorithm method)
{
    std::ofstream ofs;

    double fta; 
    double ftn; 
    double err;

    ofs.open(fname+"_"+method.method_name+".csv");
    ofs << "t,fta,ftn,err" << std::endl;
    ofs.precision(16);
    
    for(double t=1.0; t <= 10.0; t+= 1.0)
    {
        fta = ft(t);
        ftn = method.Evaluate(t, Fs);
        err = std::abs(ftn-fta)/fta;
        
        ofs << t << "," << fta << "," << ftn << "," << err << std::endl;
    }

    ofs.close();
}

void create_results_table(std::string fname)
{
    TalbotAlgorithm   talbot;
    StehfestAlgorithm stehfest;
    DeHoogAlgorithm   dehoog;

    std::vector<InverseLaplaceAlgorithm> methods{ talbot, stehfest, dehoog };

    for(InverseLaplaceAlgorithm m: methods)
    {
        output_algorithm_table(fname+"_func1", ft1, Fs1, m);
        output_algorithm_table(fname+"_func2", ft2, Fs2, m);
        output_algorithm_table(fname+"_func3", ft3, Fs3, m);
        output_algorithm_table(fname+"_func4", ft4, Fs4, m);
        output_algorithm_table(fname+"_func5", ft5, Fs5, m);
    }

}

#endif // NILT_UTILS_HEADER
