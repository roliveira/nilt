#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nilt.hpp"
#include "testfunctions.hpp"
#include "conftest.hpp"

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;


TEST_CASE("Stehfest name is really Stehfest", "[stehfest][name]")
{
    REQUIRE(std::string(nilt::Stehfest::name) == "Stehfest");
}

TEST_CASE("Stehfest throws domain_error for t <= 0", "[stehfest][domain]")
{
    nilt::Stehfest algo;
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<double>,  0.0), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<double>, -1.0), std::domain_error);
}

TEST_CASE("Stehfest throws invalid_argument for odd N or N < 2 or N > 20",
          "[stehfest][invalid]")
{
    nilt::Stehfest algo;

    for (int N : {1, 3, 5, 7, 9, 11, 13, 15, 17, 19})
    {
        CAPTURE(N);
        algo.N = N;
        REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<double>, 1.0), std::invalid_argument);
    }

    for (int N : {0, -2, -4, -6})
    {
        CAPTURE(N);
        algo.N = N;
        REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<double>, 1.0), std::invalid_argument);
    }

    for (int N : {22, 24, 26})
    {
        CAPTURE(N);
        algo.N = N;
        REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<double>, 1.0), std::invalid_argument);
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
        REQUIRE_THAT(coeff[i], Catch::Matchers::WithinRel(expected[i], TOL["STEHFEST_REL"]));
    }
}

TEST_CASE("Stehfest accepts lambda returning real", "[stehfest][callable]")
{
    nilt::Stehfest algo;
    auto Fs = [](double s) { return Fs4<double>(s); };
    double result = nilt::invert(algo, Fs, 1.0);
    REQUIRE_THAT(result, WithinRel(ft4<double>(1.0), TOL["STEHFEST_REL_TOL_SMALL"]));
}

TEST_CASE("Stehfest direct call matches free function",
          "[stehfest][api]")
{
    nilt::Stehfest algo;
    double via_free = nilt::invert(algo, Fs4<double>, 3.0);
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
        
        REQUIRE_THAT(sum, WithinAbs(0.0, TOL["STEHFEST_REL"]));

        for (double c : coeff) {
            REQUIRE(std::isfinite(c));
        }
    }
}
