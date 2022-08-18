
#ifndef NILT_PIESSENS_HEADER
#define NILT_PIESSENS_HEADER

#include <cmath>
#include <vector>

#include "invlap.hpp"


class TalbotAlgorithm: public InverseLaplaceAlgorithm
{
    public:

    int    n     = 50;
    double shift = 0.0;

    TalbotAlgorithm(void);
    double Evaluate(double t, std::complex<double> (*Fs)(std::complex<double>));
};

TalbotAlgorithm::TalbotAlgorithm(void)
: InverseLaplaceAlgorithm("Talbot") {};

double TalbotAlgorithm::Evaluate(double t, std::complex<double> (*Fs)(std::complex<double>))
{
	double theta, h;
    std::complex<double> z, dz, res;
    std::complex<double> ans(0.0, 0.0); 
	    
    h = 2*M_PI/this->n;

    for (int k=0; k<=this->n; ++k) 
    {
        theta = -M_PI + (k + 0.5)*h;
        z     = this->shift + this->n/t*(0.5017*theta/std::tan(0.6407*theta) + std::complex<double>(-0.6122, 0.2645*theta)); 
        dz    = this->n/t*(-0.5017*0.6407*theta*(std::pow(1.0/std::sin(0.6407*theta), 2.0)) + 0.5017/std::tan(0.6407*theta) + std::complex<double>(0.0, 0.2645));
        ans   = ans + std::exp(z*t)*Fs(z)*dz;
    }

    return ((h/(std::complex<double>(0.0, 2.0*M_PI)))*ans).real();
}

#endif // NILT_PIESSENS_HEADER
