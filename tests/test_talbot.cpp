#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nilt.hpp"
#include "test_tolerances.hpp"

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

TEST_CASE("Talbot inverts 1/(s+1) to exp(-t)",
          "[talbot][exp_decay]")
{
    nilt::Talbot algo;

    SECTION("t = 1.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 1.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-1.0), TALBOT_REL_TOL));
    }
    SECTION("t = 5.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 5.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-5.0), TALBOT_REL_TOL_LARGE));
    }
    SECTION("t = 10.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 10.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-10.0), TALBOT_REL_TOL_LARGE));
    }
}

TEST_CASE("Talbot inverts 1/s^2 to t",
          "[talbot][ramp]")
{
    nilt::Talbot algo;

    for (double t : {1.0, 3.0, 7.0, 10.0}) {
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs_ramp, t), WithinRel(t, TALBOT_REL_TOL));
    }
}

TEST_CASE("Talbot inverts 1/(s^2+1) to sin(t)",
          "[talbot][sin]")
{
    nilt::Talbot algo;

    for (double t : {1.0, 2.0, 5.0, 9.0}) {
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs_sin, t),
                     WithinAbs(std::sin(t), TALBOT_ABS_TOL));
    }
}

TEST_CASE("Talbot inverts s/(s^2+1) to cos(t)",
          "[talbot][cos]")
{
    nilt::Talbot algo;

    for (double t : {1.0, 3.0, 6.0}) {
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs_cos, t),
                     WithinAbs(std::cos(t), TALBOT_ABS_TOL));
    }
}

TEST_CASE("Talbot inverts 1/((s+1)^2+1) to exp(-t)sin(t)",
          "[talbot][damped_sin]")
{
    nilt::Talbot algo;

    for (double t : {1.0, 4.0, 8.0}) {
        CAPTURE(t);
        double expected = std::exp(-t) * std::sin(t);
        REQUIRE_THAT(nilt::invert(algo, Fs_damped_sin, t),
                     WithinAbs(expected, TALBOT_ABS_TOL));
    }
}

TEST_CASE("Talbot default parameters N=50 SHIFT=0", "[talbot][defaults]")
{
    nilt::Talbot algo;
    REQUIRE(algo.N == 50);
    REQUIRE(algo.SHIFT == 0.0);
}

TEST_CASE("Talbot name is Talbot", "[talbot][name]")
{
    REQUIRE(std::string(nilt::Talbot::name) == "Talbot");
}

TEST_CASE("Talbot throws domain_error for t <= 0", "[talbot][domain]")
{
    nilt::Talbot algo;
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs_exp_decay, 0.0), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs_exp_decay, -1.0), std::domain_error);
}

TEST_CASE("Talbot with N=100 inverts exp_decay",
          "[talbot][parameters]")
{
    nilt::Talbot algo;
    algo.N = 100;
    double result = nilt::invert(algo, Fs_exp_decay, 2.0);
    REQUIRE_THAT(result, WithinRel(std::exp(-2.0), TALBOT_REL_TOL_LARGE));
}

TEST_CASE("Talbot accepts complex-returning lambda", "[talbot][callable]")
{
    nilt::Talbot algo;
    auto Fs = [](C s) -> C { return 1.0 / (s + 1.0); };
    double result = nilt::invert(algo, Fs, 1.0);
    REQUIRE_THAT(result, WithinRel(std::exp(-1.0), TALBOT_REL_TOL));
}

TEST_CASE("Talbot direct call matches free function",
          "[talbot][api]")
{
    nilt::Talbot algo;
    double via_free = nilt::invert(algo, Fs_exp_decay, 3.0);
    double via_call = algo(Fs_exp_decay, 3.0);
    REQUIRE(via_free == via_call);
}
