#ifndef NLTI_UTILS_HEADER
#define NLTI_UTILS_HEADER

#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "stehfest.hpp"


double ft1(double t) { return 1.0/std::sqrt(M_PI*t); }
double ft2(double t) { return -0.57722-std::log(t); } // C = 0.57722, as it is presenteds in Stehfest (1970).
double ft3(double t) { return std::pow(t, 3.0)/6.0; }
double ft4(double t) { return std::exp(-t); }
double ft5(double t) { return std::sin(std::sqrt(2.0*t)); }
double Fs1(double s) { return 1.0/std::sqrt(s); }
double Fs2(double s) { return std::log(s)/s; }
double Fs3(double s) { return std::pow(s, -4.0); }
double Fs4(double s) { return 1.0/(s+1.0); }
double Fs5(double s) { return std::sqrt(M_PI/(2.0*std::pow(s, 3.0)))*std::exp(-1.0/(2.0*s)); }


void output_stehfest_func(std::string fname, double (*ft)(double), double (*Fs)(double))
{
    StehfestAlgorithm method;
    std::ofstream ofs;

    double fta; 
    double ftn; 
    double err;

    ofs.open(fname+".csv");
    ofs << "t,fta,ftn,err" << std::endl;
    ofs.precision(12);
    
    for(double t=1.0; t <= 10.0; t+= 1.0)
    {
        fta = ft(t);
        ftn = method.Evaluate(t, Fs);
        err = std::abs(ftn-fta)/fta;
        
        ofs << t << "," << fta << "," << ftn << "," << err << std::endl;
    }

    ofs.close();
}

void create_stehfest_table(std::string fname)
{
    output_stehfest_func(fname+"func1", ft1, Fs1);
    output_stehfest_func(fname+"func2", ft2, Fs2);
    output_stehfest_func(fname+"func3", ft3, Fs3);
    output_stehfest_func(fname+"func4", ft4, Fs4);
    output_stehfest_func(fname+"func5", ft5, Fs5);
}


#endif // NLTI_UTILS_HEADER
