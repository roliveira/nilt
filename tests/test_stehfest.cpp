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
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<double>, 0.0), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<double>, -1.0), std::domain_error);
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
    for (int N : {6, 10, 14, 18}) {
        CAPTURE(N);
        auto coeff = nilt::Stehfest::coefficients(N);
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
