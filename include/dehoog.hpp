/*
    De Hoog et al., 1982, An improved method for numerical inversion of Laplace
    transforms, SIAM J. Sci. Stat. Comput.
*/
#ifndef NILT_DEHOOG_HEADER
#define NILT_DEHOOG_HEADER

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace nilt {

class DeHoog
{
public:
    static constexpr const char* name = "DeHoog";

    int    M        = 40;       // order of approximation (number of terms)
    double T_factor = 4.0;      // period factor (T = T_factor * t)
    double tol      = 1.0e-16;  // tolerance for integration limit

    // Evaluate the inverse Laplace transform at time t.
    // Fs must be callable as Fs(std::complex<double>) -> std::complex<double>.
    template<typename F>
    double operator()(F&& Fs, double t) const
    {
        if (t <= 0.0)
            throw std::domain_error("DeHoog: t must be positive");

        int twoM = 2 * M;

        double T     = T_factor * t;
        double gamma = -0.5 * std::log(tol) / T;

        // Evaluate F(s) at quadrature points
        std::vector<std::complex<double>> Fc(twoM + 1);
        Fc[0] = 0.5 * Fs(std::complex<double>(gamma, 0.0));
        for (int i = 1; i <= twoM; ++i)
            Fc[i] = Fs(std::complex<double>(gamma, i * pi / T));

        // Quotient-difference (QD) algorithm - eq. (20) of De Hoog et al. 1982
        std::vector<std::vector<std::complex<double>>> e(twoM + 1, std::vector<std::complex<double>>(M + 1, 0.0));
        std::vector<std::vector<std::complex<double>>> q(twoM,     std::vector<std::complex<double>>(M + 1, 0.0));

        for (int i = 0; i < twoM; ++i)
            q[i][1] = Fc[i + 1] / Fc[i];

        for (int r = 1; r <= M; ++r)
        {
            for (int i = 2 * (M - r); i >= 0; --i)
                e[i][r] = q[i + 1][r] - q[i][r] + e[i + 1][r - 1];

            if (r < M)
            {
                for (int i = 2 * (M - r) - 1; i >= 0; --i)
                    q[i][r + 1] = q[i + 1][r] * e[i + 1][r] / e[i][r];
            }
        }

        // Populate d vector for continued fraction
        std::vector<std::complex<double>> d(twoM + 1);
        d[0] = Fc[0];
        for (int m = 1; m <= M; ++m)
        {
            d[2 * m - 1] = -q[0][m];
            d[2 * m]     = -e[0][m];
        }

        // Evaluate continued fraction via forward recurrence - eq. (21)
        std::complex<double> z(std::cos(pi * t / T), std::sin(pi * t / T));

        std::vector<std::complex<double>> A(twoM + 2), B(twoM + 2);
        A[0] = 0.0;  A[1] = d[0];
        B[0] = 1.0;  B[1] = 1.0;

        for (int n = 2; n <= twoM + 1; ++n)
        {
            auto dz = d[n - 1] * z;
            A[n] = A[n - 1] + dz * A[n - 2];
            B[n] = B[n - 1] + dz * B[n - 2];
        }

        // Acceleration - eqs. (23)-(24)
        auto h2M = 0.5 * (1.0 + z * (d[twoM - 1] - d[twoM]));
        auto R2M = -h2M * (1.0 - std::sqrt(1.0 + z * d[twoM] / (h2M * h2M)));

        A[twoM + 1] = A[twoM] + R2M * A[twoM - 1];
        B[twoM + 1] = B[twoM] + R2M * B[twoM - 1];

        return (1.0 / T) * std::exp(gamma * t) * (A[twoM + 1] / B[twoM + 1]).real();
    }
};

} // namespace nilt

#endif // NILT_DEHOOG_HEADER
