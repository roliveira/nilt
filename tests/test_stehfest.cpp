#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nilt.hpp"
#include "testfunctions.hpp"

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

static constexpr double STEHFEST_COEFF_TOL = 1e-2;


TEST_CASE("Stehfest name is really Stehfest", "[stehfest][name]")
{
    REQUIRE(std::string(nilt::Stehfest::name) == "Stehfest");
}

TEST_CASE("Stehfest throws domain_error for t <= 0", "[stehfest][domain]")
{
    nilt::Stehfest algo;
    REQUIRE_THROWS_AS(nilt::invert(Fs4<double>, 0.0, algo), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(Fs4<double>, -1.0, algo), std::domain_error);
}

TEST_CASE("Stehfest throws invalid_argument for odd N or N < 2 or N > 20",
          "[stehfest][invalid]")
{
    nilt::Stehfest algo;

    for (int N : {1, 3, 5, 7, 9, 11, 13, 15, 17, 19})
    {
        CAPTURE(N);
        algo.N = N;
        REQUIRE_THROWS_AS(nilt::invert(Fs4<double>, 1.0, algo), std::invalid_argument);
    }

    for (int N : {0, -2, -4, -6})
    {
        CAPTURE(N);
        algo.N = N;
        REQUIRE_THROWS_AS(nilt::invert(Fs4<double>, 1.0, algo), std::invalid_argument);
    }

    for (int N : {22, 24, 26})
    {
        CAPTURE(N);
        algo.N = N;
        REQUIRE_THROWS_AS(nilt::invert(Fs4<double>, 1.0, algo), std::invalid_argument);
    }
}

TEST_CASE("Stehfest coefficients are correct for N=18", "[stehfest][coefficients]")
{
    nilt::Stehfest algo;
    algo.N = 18;
    
    auto coeff = algo.get_coefficients();
    CAPTURE(coeff);

    REQUIRE(coeff.size() == 18);

    // Hardcoded reference coefficients for N=18 from the literature
    const std::vector<double> expected = {
         4.960317460317460e-05, -6.095734126984128e-01,  2.745940476190476e+02,
        -2.630695674603174e+04,  9.572572013888889e+05, -1.735869484583333e+07,
         1.824212226472222e+08, -1.218533288309127e+09,  5.491680025283035e+09,
        -1.736213111520684e+10,  3.945509690352738e+10, -6.526651698517500e+10,
         7.873006832822083e+10, -6.855644419612083e+10,  4.198434347505357e+10,
        -1.716093471183929e+10,  4.204550039102679e+09, -4.671722265669643e+08
    };

    for (int i = 0; i < 18; ++i)
    {
        // Use the tolerance map or a direct small value like 1e-12
        REQUIRE_THAT(coeff[i], Catch::Matchers::WithinRel(expected[i], STEHFEST_COEFF_TOL));
    }
}

TEST_CASE("Stehfest accepts lambda returning real", "[stehfest][callable]")
{
    nilt::Stehfest algo;
    auto Fs = [](double s) { return Fs4<double>(s); };
    double via_lambda = nilt::invert(Fs, 1.0, algo);
    double via_fptr   = nilt::invert(Fs4<double>, 1.0, algo);
    REQUIRE(std::isfinite(via_lambda));
    REQUIRE_THAT(via_lambda, WithinRel(via_fptr, 1e-12));
}

TEST_CASE("Stehfest direct call matches free function",
          "[stehfest][api]")
{
    nilt::Stehfest algo;
    double via_free = nilt::invert(Fs4<double>, 3.0, algo);
    double via_call = algo(Fs4<double>, 3.0);
    REQUIRE(via_free == via_call);
}

TEST_CASE("Stehfest coefficients sum to zero for even N",
          "[stehfest][coefficients]")
{
    for (int N : {2, 4, 6, 8, 10, 12, 14, 16, 18, 20}) {
        CAPTURE(N);
        
        nilt::Stehfest algo;
        algo.N = N;
        auto coeff = algo.get_coefficients();

        REQUIRE(coeff.size() == static_cast<std::size_t>(N));
        
        double sum = 0.0;
        for (double c : coeff)
            sum += c;
        
        REQUIRE_THAT(sum, WithinAbs(0.0, STEHFEST_COEFF_TOL));

        for (double c : coeff) {
            REQUIRE(std::isfinite(c));
        }
    }
}

TEST_CASE("Stehfest array eval_batched matches operator() per t",
          "[stehfest][eval_batched]")
{
    nilt::Stehfest algo;
    const std::vector<double> t = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
    std::vector<double> out(t.size());

    algo.eval_batched(
        [](const double* s, double* f, int n) {
            for (int i = 0; i < n; ++i) f[i] = Fs4<double>(s[i]);
        },
        t.data(), out.data(), static_cast<int>(t.size()));

    for (size_t i = 0; i < t.size(); ++i) {
        CAPTURE(t[i]);
        REQUIRE_THAT(out[i], WithinRel(algo(Fs4<double>, t[i]), 1e-12));
    }
}

TEST_CASE("Stehfest array eval_batched invokes Fs_batched exactly once",
          "[stehfest][eval_batched]")
{
    nilt::Stehfest algo;
    const std::vector<double> t = {0.5, 1.0, 2.0, 3.0};
    std::vector<double> out(t.size());
    int call_count = 0;
    int total_points = 0;

    algo.eval_batched(
        [&](const double* s, double* f, int n) {
            ++call_count;
            total_points = n;
            for (int i = 0; i < n; ++i) f[i] = Fs4<double>(s[i]);
        },
        t.data(), out.data(), static_cast<int>(t.size()));

    REQUIRE(call_count == 1);
    REQUIRE(total_points == static_cast<int>(t.size()) * algo.N);
}

TEST_CASE("Stehfest array eval_batched rejects non-positive t",
          "[stehfest][eval_batched][domain]")
{
    nilt::Stehfest algo;
    const std::vector<double> t_bad = {1.0, 2.0, 0.0, 3.0};
    std::vector<double> out(t_bad.size());
    auto Fb = [](const double* s, double* f, int n) {
        for (int i = 0; i < n; ++i) f[i] = Fs4<double>(s[i]);
    };
    REQUIRE_THROWS_AS(
        algo.eval_batched(Fb, t_bad.data(), out.data(), 4),
        std::domain_error);
}
