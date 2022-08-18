
#ifndef NILT_TESTFUNCTIONS_HEADER
#define NILT_TESTFUNCTIONS_HEADER

#include <cmath>


template<typename T> T ft1(T t) { return 1.0/std::sqrt(M_PI*t); }
template<typename T> T ft2(T t) { return -0.57722-std::log(t); } // C = 0.57722, as it is presented in Stehfest (1970).
template<typename T> T ft3(T t) { return std::pow(t, 3.0)/6.0; }
template<typename T> T ft4(T t) { return std::exp(-t); }
template<typename T> T ft5(T t) { return std::sin(std::sqrt(2.0*t)); }

template<typename T> T Fs1(T s) { return 1.0/std::sqrt(s); }
template<typename T> T Fs2(T s) { return std::log(s)/s; }
template<typename T> T Fs3(T s) { return std::pow(s, -4.0); }
template<typename T> T Fs4(T s) { return 1.0/(s+1.0); }
template<typename T> T Fs5(T s) { return std::sqrt(M_PI/(2.0*std::pow(s, 3.0)))*std::exp(-1.0/(2.0*s)); }


#endif // NILT_TESTFUNCTIONS_HEADER
