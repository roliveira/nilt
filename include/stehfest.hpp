/*
    Harald Stehfest. 1970. Algorithm 368: Numerical inversion of Laplace transforms [D5]. 
    Commun. ACM 13, 1 (Jan. 1970), 47-49. DOI:https://doi.org/10.1145/361953.361969
*/
#ifndef NILT_STEHFEST_HEADER
#define NILT_STEHFEST_HEADER

#include <cmath>
#include <stdexcept>
#include <vector>

namespace nilt {

class Stehfest
{
public:
    static constexpr const char* name = "Stehfest";

    int N = 18;  // number of terms (must be even)

    // Evaluate the inverse Laplace transform at time t.
    // Fs must be callable as Fs(double) -> double.
    template<typename F>
    double operator()(F&& Fs, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("Stehfest: t must be positive");

        auto coeff = coefficients(N);
        double ln2t = std::log(2.0) / t;
        double s = 0.0;
        double y = 0.0;

        for (double c : coeff)
        {
            s += ln2t;
            y += c * Fs(s);
        }

        return ln2t * y;
    }

    // Compute Stehfest coefficients V_i for i=1..N.
    static std::vector<double> coefficients(int N)
    {
        int N2 = N / 2;
        std::vector<double> V(N);

        int sign = (N2 % 2 != 0) ? -1 : 1;

        for (int i = 0; i < N; ++i)
        {
            int kmin = (i + 2) / 2;
            int kmax = std::min(i + 1, N2);

            double sum = 0.0;
            sign = -sign;

            for (int k = kmin; k <= kmax; ++k)
            {
                sum += std::pow(static_cast<double>(k), N2)
                     * factorial(2 * k)
                     / (factorial(k) * factorial(k - 1)
                        * factorial(N2 - k) * factorial(i + 1 - k)
                        * factorial(2 * k - i - 1));
            }

            V[i] = sign * sum;
        }

        return V;
    }

private:
    static double factorial(int n)
    {
        double x = 1.0;
        for (int i = 2; i <= n; ++i)
            x *= i;
        return x;
    }
};

} // namespace nilt

#endif // NILT_STEHFEST_HEADER
