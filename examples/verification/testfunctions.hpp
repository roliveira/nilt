
#ifndef NILT_TESTFUNCTIONS_HEADER
#define NILT_TESTFUNCTIONS_HEADER

#include <cmath>

#include "nilt.hpp"


// Stehfest (1970) test functions

// f1: f(t) = 1/sqrt(pi*t),        F(s) = 1/sqrt(s)
template<typename T> T ft1(T t) { return 1.0/std::sqrt(nilt::pi*t); }
template<typename T> T Fs1(T s) { return 1.0/std::sqrt(s); }

// f2: f(t) = -gamma - ln(t),       F(s) = ln(s)/s
template<typename T> T ft2(T t) { return -0.57722-std::log(t); }
template<typename T> T Fs2(T s) { return std::log(s)/s; }

// f3: f(t) = t^3/6,                F(s) = 1/s^4
template<typename T> T ft3(T t) { return std::pow(t, 3.0)/6.0; }
template<typename T> T Fs3(T s) { return std::pow(s, -4.0); }

// f4: f(t) = exp(-t),              F(s) = 1/(s+1)
template<typename T> T ft4(T t) { return std::exp(-t); }
template<typename T> T Fs4(T s) { return 1.0/(s+1.0); }

// f5: f(t) = sin(sqrt(2t)),        F(s) = sqrt(pi/(2s^3)) * exp(-1/(2s))
template<typename T> T ft5(T t) { return std::sin(std::sqrt(2.0*t)); }
template<typename T> T Fs5(T s) { return std::sqrt(nilt::pi/(2.0*std::pow(s, 3.0)))*std::exp(-1.0/(2.0*s)); }

// Abate & Whitt (2006) test functions

// f6: f(t) = t,                    F(s) = 1/s^2
template<typename T> T ft6(T t) { return t; }
template<typename T> T Fs6(T s) { return 1.0/(s*s); }

// f7: f(t) = t*exp(-t),            F(s) = 1/(s+1)^2
template<typename T> T ft7(T t) { return t*std::exp(-t); }
template<typename T> T Fs7(T s) { return 1.0/((s+1.0)*(s+1.0)); }

// f8: f(t) = sin(t),               F(s) = 1/(s^2+1)
template<typename T> T ft8(T t) { return std::sin(t); }
template<typename T> T Fs8(T s) { return 1.0/(s*s+1.0); }

// f9: f(t) = cos(t),               F(s) = s/(s^2+1)
template<typename T> T ft9(T t) { return std::cos(t); }
template<typename T> T Fs9(T s) { return s/(s*s+1.0); }

// f10: f(t) = exp(-t)*sin(t),      F(s) = 1/((s+1)^2+1)
template<typename T> T ft10(T t) { return std::exp(-t)*std::sin(t); }
template<typename T> T Fs10(T s) { return 1.0/((s+1.0)*(s+1.0)+1.0); }


#endif // NILT_TESTFUNCTIONS_HEADER
