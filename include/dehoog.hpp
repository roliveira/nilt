/*
    De Hoog et al., 1982, An improved method for numerical inversion of Laplace
    transforms, SIAM J. Sci. Stat. Comput.

    Laplace Inversion Source Code - Copyright © 2010 James R. Craig, University of Waterloo
*/
#ifndef NILT_DEHOOG_HEADER
#define NILT_DEHOOG_HEADER


#include <cmath>
#include <complex>

#include "invlap.hpp"


class DeHoogAlgorithm: public InverseLaplaceAlgorithm
{
    public:

    int    MAX_LAPORDER = 60;
    int    M            = 40;       // order of Taylor Expansion (must be less than MAX_LAPORDER)
    double DeHoogFactor = 4.0;      // DeHoog time factor
    double tol          = 1.0e-16;  // tolerance
    
    DeHoogAlgorithm(void);

    double Evaluate(double t, std::complex<double> (*Fs)(std::complex<double>));
};

DeHoogAlgorithm::DeHoogAlgorithm(void)
: InverseLaplaceAlgorithm("DeHoog") {};

double DeHoogAlgorithm::Evaluate(double t, std::complex<double> (*Fs)(std::complex<double>))
{
    int    i, n, m, r;          // counters and intermediate array indices
    double T;                   // Period of DeHoog Inversion formula
    double gamma;               // Integration limit parameter
    
    std::complex<double> h2M, R2M, z, dz, s;  // Temporary variables

    std::complex<double> Fctrl [2*this->MAX_LAPORDER+1];
    std::complex<double> e     [2*this->MAX_LAPORDER  ][this->MAX_LAPORDER];
    std::complex<double> q     [2*this->MAX_LAPORDER  ][this->MAX_LAPORDER];
    std::complex<double> d     [2*this->MAX_LAPORDER+1];
    std::complex<double> A     [2*this->MAX_LAPORDER+2]; 
    std::complex<double> B     [2*this->MAX_LAPORDER+2]; 

    T     = this->DeHoogFactor*t;        // calculate the period
    gamma = -0.5*std::log(this->tol)/T;  // and the integration limits

    // Calculate F(s) at evaluation points gamma+IM*i*PI/T for i=0 to 2*M-1
    // This is likely the most time consuming portion of the DeHoog algorithm
    Fctrl[0] = 0.5*Fs(gamma); 

    for (i=1; i<=2*this->M; ++i)
    {
        s = std::complex<double>(gamma, i*M_PI/T);
        Fctrl[i] = Fs(s);
    }

    // Evaluate e and q 
    // eq. (20) of de Hoog et al 1982
    for (i=0; i<2*this->M; ++i)
    {
        e[i][0] = 0.0;
        q[i][1] = Fctrl[i+1]/Fctrl[i];
    }

    e[2*this->M][0] = 0.0;

    for (r=1; r<=this->M-1; r++) // one minor correction - does not work for r<=M, as suggested in paper
    {
        for (i=2*(this->M-r); i>=0; --i)
        {     
            if ((i<2*(this->M-r)) && (r>1)) 
            {
                q[i][r] = q[i+1][r-1]*e[i+1][r-1]/e[i][r-1];
            }

            e[i][r]=q[i+1][r]-q[i][r]+e[i+1][r-1];
        }
    }

    // Populate d vector
    d[0] = Fctrl[0];

    for (m=1; m<=this->M; ++m)
    {
        d[2*m-1] = -q[0][m];
        d[2*m  ] = -e[0][m];
    }

    // Evaluate A, B
    // eq. (21) in De Hoog et al.
    z = std::complex<double>(std::cos(M_PI*t/T), std::sin(M_PI*t/T));

    A[0] = 0.0;   // A_{-1}, B_{-1} in De Hoog
    A[1] = d[0]; 
    B[0] = 1.0; 
    B[1] = 1.0;

    for (n=2; n<=2*this->M+1; ++n) 
    {
        dz   = d[n-1]*z; 
        A[n] = A[n-1]+dz*A[n-2];
        B[n] = B[n-1]+dz*B[n-2];
    }

    // eq. (23) in De Hoog et al.
    h2M = 0.5*(1.0+z*(d[2*this->M-1]-d[2*this->M]));
    R2M = -h2M*(1.0-std::sqrt(1.0+(z*d[2*this->M]/h2M/h2M)));

    // eq. (24) in De Hoog et al.
    A[2*this->M+1]=A[2*this->M]+R2M*A[2*this->M-1];
    B[2*this->M+1]=B[2*this->M]+R2M*B[2*this->M-1];

    // Final result: A[2*M]/B[2*M] = sum [F(gamma+itheta)*exp(itheta)]
    return 1.0/T*std::exp(gamma*t)*(A[2*this->M+1]/B[2*this->M+1]).real();
}

#endif // NILT_DEHOOG_HEADER
