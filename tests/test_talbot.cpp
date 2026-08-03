#include <complex>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nilt.hpp"
#include "testfunctions.hpp"


TEST_CASE("Talbot default parameters N=50 SHIFT=0", "[talbot][defaults]")
{
    nilt::Talbot algo;
    REQUIRE(algo.N == 50);
    REQUIRE(algo.SHIFT == 0.0);
}

TEST_CASE("Talbot name is really Talbot", "[talbot][name]")
{
    REQUIRE(std::string(nilt::Talbot::name) == "Talbot");
}

TEST_CASE("Talbot throws domain_error for t <= 0", "[talbot][domain]")
{
    nilt::Talbot algo;
    REQUIRE_THROWS_AS(nilt::invert(Fs4<std::complex<double>>, 0.0, algo), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(Fs4<std::complex<double>>, -1.0, algo), std::domain_error);
}

TEST_CASE("Talbot invert with N=64 (table boundary)", "[talbot][parameters]")
{
    nilt::Talbot algo;
    algo.N = 64;
    double result = nilt::invert(Fs4<std::complex<double>>, 2.0, algo);
    REQUIRE(std::isfinite(result));
}

TEST_CASE("Talbot invert with N=100 (runtime fallback)", "[talbot][parameters]")
{
    nilt::Talbot algo;
    algo.N = 100;
    double result = nilt::invert(Fs4<std::complex<double>>, 2.0, algo);
    REQUIRE(std::isfinite(result));
}

TEST_CASE("Talbot accepts complex-returning lambda", "[talbot][callable]")
{
    nilt::Talbot algo;
    auto Fs = [](std::complex<double> s) -> std::complex<double> { return Fs4<std::complex<double>>(s); };
    double via_lambda = nilt::invert(Fs, 1.0, algo);
    double via_fptr   = nilt::invert(Fs4<std::complex<double>>, 1.0, algo);
    REQUIRE(std::isfinite(via_lambda));
    REQUIRE_THAT(via_lambda, Catch::Matchers::WithinRel(via_fptr, 1e-12));
}

TEST_CASE("Talbot direct call matches free function", "[talbot][api]")
{
    nilt::Talbot algo;
    double via_free = nilt::invert(Fs4<std::complex<double>>, 3.0, algo);
    double via_call = algo(Fs4<std::complex<double>>, 3.0);
    REQUIRE(via_free == via_call);
}

TEST_CASE("Talbot array eval_batched matches operator() per t",
          "[talbot][eval_batched]")
{
    using C = std::complex<double>;
    nilt::Talbot algo;
    const std::vector<double> t = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
    std::vector<double> out(t.size());

    algo.eval_batched(
        [](const C* s, C* f, int n) {
            for (int i = 0; i < n; ++i) f[i] = Fs4<C>(s[i]);
        },
        t.data(), out.data(), static_cast<int>(t.size()));

    for (size_t i = 0; i < t.size(); ++i) {
        CAPTURE(t[i]);
        double ref = algo(Fs4<C>, t[i]);
        REQUIRE(std::abs(out[i] - ref) < 1e-12 * std::max(std::abs(ref), 1.0));
    }
}

TEST_CASE("Talbot array eval_batched invokes Fs_batched exactly once",
          "[talbot][eval_batched]")
{
    using C = std::complex<double>;
    nilt::Talbot algo;
    const std::vector<double> t = {0.5, 1.0, 2.0, 3.0};
    std::vector<double> out(t.size());
    int call_count = 0;
    int total_points = 0;

    algo.eval_batched(
        [&](const C* s, C* f, int n) {
            ++call_count;
            total_points = n;
            for (int i = 0; i < n; ++i) f[i] = Fs4<C>(s[i]);
        },
        t.data(), out.data(), static_cast<int>(t.size()));

    REQUIRE(call_count == 1);
    REQUIRE(total_points == static_cast<int>(t.size()) * algo.N);
}

TEST_CASE("Talbot array eval_batched works for runtime-fallback N=100",
          "[talbot][eval_batched]")
{
    using C = std::complex<double>;
    nilt::Talbot algo;
    algo.N = 100;
    const std::vector<double> t = {0.5, 1.0, 2.0};
    std::vector<double> out(t.size());

    algo.eval_batched(
        [](const C* s, C* f, int n) {
            for (int i = 0; i < n; ++i) f[i] = Fs4<C>(s[i]);
        },
        t.data(), out.data(), static_cast<int>(t.size()));

    for (size_t i = 0; i < t.size(); ++i) {
        CAPTURE(t[i]);
        double ref = algo(Fs4<C>, t[i]);
        REQUIRE(std::abs(out[i] - ref) < 1e-10 * std::max(std::abs(ref), 1.0));
    }
}

TEST_CASE("Talbot array eval_batched rejects non-positive t",
          "[talbot][eval_batched][domain]")
{
    using C = std::complex<double>;
    nilt::Talbot algo;
    const std::vector<double> t_bad = {1.0, -0.5, 2.0};
    std::vector<double> out(t_bad.size());
    auto Fb = [](const C* s, C* f, int n) {
        for (int i = 0; i < n; ++i) f[i] = Fs4<C>(s[i]);
    };
    REQUIRE_THROWS_AS(
        algo.eval_batched(Fb, t_bad.data(), out.data(), 3),
        std::domain_error);
}

TEST_CASE("Talbot yields finite output at every N in [1, 128]",
          "[talbot][nan][singularity]")
{
    using C = std::complex<double>;
    // Sweep every N in the table range plus the runtime-fallback boundaries.
    // Odd N values (e.g. 7, 9, 65) hit theta = 0 at the midpoint k = N/2 and
    // trip the removable-singularity branch; the test locks in that both the
    // table path and the fallback path return finite values there.
    nilt::Talbot algo;
    const double t = 1.5;
    const double ref = std::exp(-t);  // Fs4 has f(t) = exp(-t)
    for (int N = 1; N <= 128; ++N) {
        CAPTURE(N);
        algo.N = N;
        double r = algo(Fs4<C>, t);
        REQUIRE(std::isfinite(r));
        // Coarse sanity: for N >= 8 the Talbot method should be well within
        // 1% of the analytic value.  Below 8 the algorithm is unusable but
        // must at least not NaN.
        if (N >= 8 && N <= 64) {
            REQUIRE(std::abs(r - ref) < 1e-2);
        }
    }
}
