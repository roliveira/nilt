#include <complex>

#include <catch2/catch_test_macros.hpp>

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
    REQUIRE(via_lambda == via_fptr);
}

TEST_CASE("DeHoog direct call matches free function",
          "[dehoog][api]")
{
    nilt::DeHoog algo;
    double via_free = nilt::invert(algo, Fs4<std::complex<double>>, 3.0);
    double via_call = algo(Fs4<std::complex<double>>, 3.0);
    REQUIRE(via_free == via_call);
}
