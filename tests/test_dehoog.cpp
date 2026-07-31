#include <complex>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nilt.hpp"
#include "testfunctions.hpp"



TEST_CASE("DeHoog default parameters M=40 T_FACTOR=4.0 TOL=1e-16",
          "[dehoog][defaults]")
{
    nilt::DeHoog algo;
    REQUIRE(algo.M == 40);
    REQUIRE(algo.T_FACTOR == 4.0);
    REQUIRE(algo.TOL == 1e-16);
}

TEST_CASE("DeHoog name is DeHoog", "[dehoog][name]")
{
    REQUIRE(std::string(nilt::DeHoog::name) == "DeHoog");
}

TEST_CASE("DeHoog throws domain_error for t <= 0", "[dehoog][domain]")
{
    nilt::DeHoog algo;
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<std::complex<double>>, 0.0), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<std::complex<double>>, -1.0), std::domain_error);
}

TEST_CASE("DeHoog with M=60 produces finite result",
          "[dehoog][parameters]")
{
    nilt::DeHoog algo;
    algo.M = 60;
    double result = nilt::invert(algo, Fs4<std::complex<double>>, 2.0);
    REQUIRE(std::isfinite(result));
}

TEST_CASE("DeHoog accepts complex-returning lambda", "[dehoog][callable]")
{
    nilt::DeHoog algo;
    auto Fs = [](std::complex<double> s) -> std::complex<double> { return 1.0 / (s + 1.0); };
    double via_lambda = nilt::invert(algo, Fs, 1.0);
    double via_fptr   = nilt::invert(algo, Fs4<std::complex<double>>, 1.0);
    REQUIRE(std::isfinite(via_lambda));
    REQUIRE_THAT(via_lambda, Catch::Matchers::WithinRel(via_fptr, 1e-12));
}

TEST_CASE("DeHoog direct call matches free function",
          "[dehoog][api]")
{
    nilt::DeHoog algo;
    double via_free = nilt::invert(algo, Fs4<std::complex<double>>, 3.0);
    double via_call = algo(Fs4<std::complex<double>>, 3.0);
    REQUIRE(via_free == via_call);
}

TEST_CASE("DeHoog array eval_batched matches operator() per t",
          "[dehoog][eval_batched]")
{
    using C = std::complex<double>;
    nilt::DeHoog algo;
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
        REQUIRE(std::abs(out[i] - ref) < 1e-10 * std::max(std::abs(ref), 1.0));
    }
}

TEST_CASE("DeHoog array eval_batched invokes Fs_batched exactly once",
          "[dehoog][eval_batched]")
{
    using C = std::complex<double>;
    nilt::DeHoog algo;
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
    REQUIRE(total_points == static_cast<int>(t.size()) * (2 * algo.M + 1));
}

TEST_CASE("DeHoog array eval_batched rejects non-positive t",
          "[dehoog][eval_batched][domain]")
{
    using C = std::complex<double>;
    nilt::DeHoog algo;
    const std::vector<double> t_bad = {1.0, 2.0, 0.0};
    std::vector<double> out(t_bad.size());
    auto Fb = [](const C* s, C* f, int n) {
        for (int i = 0; i < n; ++i) f[i] = Fs4<C>(s[i]);
    };
    REQUIRE_THROWS_AS(
        algo.eval_batched(Fb, t_bad.data(), out.data(), 3),
        std::domain_error);
}
