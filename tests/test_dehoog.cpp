#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nilt.hpp"

#include <complex>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;
using C = std::complex<double>;

// F(s) = 1/(s+1)  =>  f(t) = exp(-t)
static C Fs_exp_decay(C s) { return 1.0 / (s + 1.0); }

// F(s) = 1/s^2  =>  f(t) = t
static C Fs_ramp(C s) { return 1.0 / (s * s); }

// F(s) = 1/(s^2+1)  =>  f(t) = sin(t)
static C Fs_sin(C s) { return 1.0 / (s * s + 1.0); }

// F(s) = s/(s^2+1)  =>  f(t) = cos(t)
static C Fs_cos(C s) { return s / (s * s + 1.0); }

// F(s) = 1/((s+1)^2+1)  =>  f(t) = exp(-t)*sin(t)
static C Fs_damped_sin(C s) { return 1.0 / ((s + 1.0) * (s + 1.0) + 1.0); }

TEST_CASE("DeHoog inverts 1/(s+1) to exp(-t) within 1e-12",
          "[dehoog][exp_decay]")
{
    nilt::DeHoog algo;

    SECTION("t = 1.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 1.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-1.0), 1e-12));
    }
    SECTION("t = 5.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 5.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-5.0), 1e-12));
    }
    SECTION("t = 10.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 10.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-10.0), 1e-9));
    }
}

TEST_CASE("DeHoog inverts 1/s^2 to t within 1e-12",
          "[dehoog][ramp]")
{
    nilt::DeHoog algo;

    for (double t : {1.0, 3.0, 7.0, 10.0}) {
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs_ramp, t), WithinRel(t, 1e-12));
    }
}

TEST_CASE("DeHoog inverts 1/(s^2+1) to sin(t) within 1e-12",
          "[dehoog][sin]")
{
    nilt::DeHoog algo;

    for (double t : {1.0, 2.0, 5.0, 9.0}) {
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs_sin, t),
                     WithinAbs(std::sin(t), 1e-12));
    }
}

TEST_CASE("DeHoog inverts s/(s^2+1) to cos(t) within 1e-12",
          "[dehoog][cos]")
{
    nilt::DeHoog algo;

    for (double t : {1.0, 3.0, 6.0}) {
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs_cos, t),
                     WithinAbs(std::cos(t), 1e-12));
    }
}

TEST_CASE("DeHoog inverts 1/((s+1)^2+1) to exp(-t)sin(t) within 1e-12",
          "[dehoog][damped_sin]")
{
    nilt::DeHoog algo;

    for (double t : {1.0, 4.0, 8.0}) {
        CAPTURE(t);
        double expected = std::exp(-t) * std::sin(t);
        REQUIRE_THAT(nilt::invert(algo, Fs_damped_sin, t),
                     WithinAbs(expected, 1e-12));
    }
}

TEST_CASE("DeHoog default parameters M=40 T_factor=4.0 tol=1e-16",
          "[dehoog][defaults]")
{
    nilt::DeHoog algo;
    REQUIRE(algo.M == 40);
    REQUIRE(algo.T_factor == 4.0);
    REQUIRE(algo.tol == 1e-16);
}

TEST_CASE("DeHoog name is DeHoog", "[dehoog][name]")
{
    REQUIRE(std::string(nilt::DeHoog::name) == "DeHoog");
}

TEST_CASE("DeHoog throws domain_error for t <= 0", "[dehoog][domain]")
{
    nilt::DeHoog algo;
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs_exp_decay, 0.0), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs_exp_decay, -1.0), std::domain_error);
}

TEST_CASE("DeHoog with M=60 inverts exp_decay within 1e-11",
          "[dehoog][parameters]")
{
    nilt::DeHoog algo;
    algo.M = 60;
    double result = nilt::invert(algo, Fs_exp_decay, 2.0);
    REQUIRE_THAT(result, WithinRel(std::exp(-2.0), 1e-11));
}

TEST_CASE("DeHoog accepts complex-returning lambda", "[dehoog][callable]")
{
    nilt::DeHoog algo;
    auto Fs = [](C s) -> C { return 1.0 / (s + 1.0); };
    double result = nilt::invert(algo, Fs, 1.0);
    REQUIRE_THAT(result, WithinRel(std::exp(-1.0), 1e-12));
}

TEST_CASE("DeHoog direct call matches free function",
          "[dehoog][api]")
{
    nilt::DeHoog algo;
    double via_free = nilt::invert(algo, Fs_exp_decay, 3.0);
    double via_call = algo(Fs_exp_decay, 3.0);
    REQUIRE(via_free == via_call);
}
