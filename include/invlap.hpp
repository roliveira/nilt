
#ifndef NILT_INVLAP_HEADER
#define NILT_INVLAP_HEADER


#include <complex>
#include <string>


class InverseLaplaceAlgorithm
{
    public:

    std::string method_name;

    InverseLaplaceAlgorithm(std::string);

    virtual double Evaluate(double t, double (*Fs)(double)) { return 0.0; };
    virtual double Evaluate(double t, std::complex<double> (*Fs)(std::complex<double>)) { return 0.0; };
};

InverseLaplaceAlgorithm::InverseLaplaceAlgorithm(std::string name)
: method_name(name) {};


#endif // NILT_INVLAP_HEADER
