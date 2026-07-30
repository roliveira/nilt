#include <complex>

#include <catch2/catch_test_macros.hpp>

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
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<std::complex<double>>, 0.0), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs4<std::complex<double>>, -1.0), std::domain_error);
}

TEST_CASE("Talbot invert with N=64 (table boundary)", "[talbot][parameters]")
{
    nilt::Talbot algo;
    algo.N = 64;
    double result = nilt::invert(algo, Fs4<std::complex<double>>, 2.0);
    REQUIRE(std::isfinite(result));
}

TEST_CASE("Talbot invert with N=100 (runtime fallback)", "[talbot][parameters]")
{
    nilt::Talbot algo;
    algo.N = 100;
    double result = nilt::invert(algo, Fs4<std::complex<double>>, 2.0);
    REQUIRE(std::isfinite(result));
}

TEST_CASE("Talbot accepts complex-returning lambda", "[talbot][callable]")
{
    nilt::Talbot algo;
    auto Fs = [](std::complex<double> s) -> std::complex<double> { return Fs4<std::complex<double>>(s); };
    double via_lambda = nilt::invert(algo, Fs, 1.0);
    double via_fptr   = nilt::invert(algo, Fs4<std::complex<double>>, 1.0);
    REQUIRE(std::isfinite(via_lambda));
    REQUIRE(via_lambda == via_fptr);
}

TEST_CASE("Talbot direct call matches free function", "[talbot][api]")
{
    nilt::Talbot algo;
    double via_free = nilt::invert(algo, Fs4<std::complex<double>>, 3.0);
    double via_call = algo(Fs4<std::complex<double>>, 3.0);
    REQUIRE(via_free == via_call);
}
