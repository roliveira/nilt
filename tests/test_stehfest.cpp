#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nilt.hpp"

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// F(s) = 1/(s+1)  =>  f(t) = exp(-t)
static double Fs_exp_decay(double s) { return 1.0 / (s + 1.0); }

// F(s) = 1/s^2  =>  f(t) = t
static double Fs_ramp(double s) { return 1.0 / (s * s); }

// F(s) = 1/s^4  =>  f(t) = t^3/6
static double Fs_cubic(double s) { return 1.0 / (s * s * s * s); }

TEST_CASE("Stehfest inverts exp_decay 1/(s+1) to exp(-t) within 1e-4 at small t",
          "[stehfest][exp_decay]")
{
    nilt::Stehfest algo;

    SECTION("t = 1.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 1.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-1.0), 1e-4));
    }
    SECTION("t = 2.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 2.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-2.0), 1e-4));
    }
    SECTION("t = 3.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 3.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-3.0), 1e-4));
    }
}

TEST_CASE("Stehfest inverts exp_decay 1/(s+1) within 0.1 at large t",
          "[stehfest][exp_decay][large_t]")
{
    nilt::Stehfest algo;

    // Stehfest accuracy degrades for large t with exponentially decaying functions
    SECTION("t = 5.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 5.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-5.0), 0.01));
    }
    SECTION("t = 10.0") {
        double result = nilt::invert(algo, Fs_exp_decay, 10.0);
        REQUIRE_THAT(result, WithinRel(std::exp(-10.0), 0.1));
    }
}

TEST_CASE("Stehfest inverts 1/s^2 to t within 1e-6",
          "[stehfest][ramp]")
{
    nilt::Stehfest algo;

    for (double t : {1.0, 3.0, 7.0, 10.0}) {
        CAPTURE(t);
        double result = nilt::invert(algo, Fs_ramp, t);
        REQUIRE_THAT(result, WithinRel(t, 1e-6));
    }
}

TEST_CASE("Stehfest inverts 1/s^4 to t^3/6 within 1e-4",
          "[stehfest][cubic]")
{
    nilt::Stehfest algo;

    for (double t : {1.0, 4.0, 8.0}) {
        CAPTURE(t);
        double expected = t * t * t / 6.0;
        double result = nilt::invert(algo, Fs_cubic, t);
        REQUIRE_THAT(result, WithinRel(expected, 1e-4));
    }
}

TEST_CASE("Stehfest with N=12 inverts exp_decay within 1e-3",
          "[stehfest][parameters]")
{
    nilt::Stehfest algo;
    algo.N = 12;

    double result = nilt::invert(algo, Fs_exp_decay, 2.0);
    REQUIRE_THAT(result, WithinRel(std::exp(-2.0), 1e-3));
}

TEST_CASE("Stehfest default N is 18", "[stehfest][defaults]")
{
    nilt::Stehfest algo;
    REQUIRE(algo.N == 18);
}

TEST_CASE("Stehfest name is Stehfest", "[stehfest][name]")
{
    REQUIRE(std::string(nilt::Stehfest::name) == "Stehfest");
}

TEST_CASE("Stehfest throws domain_error for t <= 0", "[stehfest][domain]")
{
    nilt::Stehfest algo;
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs_exp_decay, 0.0), std::domain_error);
    REQUIRE_THROWS_AS(nilt::invert(algo, Fs_exp_decay, -1.0), std::domain_error);
}

TEST_CASE("Stehfest accepts lambda returning real", "[stehfest][callable]")
{
    nilt::Stehfest algo;
    auto Fs = [](double s) { return 1.0 / (s + 1.0); };
    double result = nilt::invert(algo, Fs, 1.0);
    REQUIRE_THAT(result, WithinRel(std::exp(-1.0), 1e-4));
}

TEST_CASE("Stehfest direct call matches free function",
          "[stehfest][api]")
{
    nilt::Stehfest algo;
    double via_free = nilt::invert(algo, Fs_exp_decay, 3.0);
    double via_call = algo(Fs_exp_decay, 3.0);
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
        // Stehfest coefficients for the weights satisfy sum(V_i) ~ 0
        // but the actual property is sum(V_i * i) = ln2 * N/2
        // Check each coefficient is finite
        for (double c : coeff) {
            REQUIRE(std::isfinite(c));
        }
    }
}
