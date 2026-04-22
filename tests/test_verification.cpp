#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nilt.hpp"

#include <cmath>
#include <complex>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;
using C = std::complex<double>;

// Cross-algorithm consistency: all three methods agree on the same F(s)

TEST_CASE("All three methods agree on 1/(s+1) -> exp(-t) at t=2",
          "[cross][exp_decay]")
{
    double expected = std::exp(-2.0);

    nilt::Stehfest stehfest;
    nilt::Talbot   talbot;
    nilt::DeHoog   dehoog;

    double s_val = nilt::invert(stehfest, [](double s) { return 1.0 / (s + 1.0); }, 2.0);
    double t_val = nilt::invert(talbot, [](C s) -> C { return 1.0 / (s + 1.0); }, 2.0);
    double d_val = nilt::invert(dehoog, [](C s) -> C { return 1.0 / (s + 1.0); }, 2.0);

    // Stehfest is less precise (1e-4), Talbot and DeHoog much better
    REQUIRE_THAT(s_val, WithinRel(expected, 1e-4));
    REQUIRE_THAT(t_val, WithinRel(expected, 1e-10));
    REQUIRE_THAT(d_val, WithinRel(expected, 1e-12));
}

// Verification suite: 10 standard test functions at t=1..10
// Tests exact numerical values against analytical solutions.

// Tolerances per method (relative error bounds observed in benchmarks)
static constexpr double TOL_STEHFEST = 0.1;
static constexpr double TOL_TALBOT   = 1e-6;
static constexpr double TOL_DEHOOG   = 1e-9;

// func1: f(t) = 1/sqrt(pi*t),  F(s) = 1/sqrt(s)

TEST_CASE("Verification func1: Stehfest inverts 1/sqrt(s) at t=1..10",
          "[verification][func1][stehfest]")
{
    nilt::Stehfest algo;
    auto Fs = [](double s) { return 1.0 / std::sqrt(s); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        double expected = 1.0 / std::sqrt(nilt::pi * t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinRel(expected, TOL_STEHFEST));
    }
}

TEST_CASE("Verification func1: Talbot inverts 1/sqrt(s) at t=1..10",
          "[verification][func1][talbot]")
{
    nilt::Talbot algo;
    auto Fs = [](C s) -> C { return 1.0 / std::sqrt(s); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        double expected = 1.0 / std::sqrt(nilt::pi * t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinRel(expected, TOL_TALBOT));
    }
}

TEST_CASE("Verification func1: DeHoog inverts 1/sqrt(s) at t=1..10",
          "[verification][func1][dehoog]")
{
    nilt::DeHoog algo;
    auto Fs = [](C s) -> C { return 1.0 / std::sqrt(s); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        double expected = 1.0 / std::sqrt(nilt::pi * t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinRel(expected, TOL_DEHOOG));
    }
}

// func4: f(t) = exp(-t),  F(s) = 1/(s+1)

TEST_CASE("Verification func4: Stehfest inverts 1/(s+1) at t=1..10",
          "[verification][func4][stehfest]")
{
    nilt::Stehfest algo;
    auto Fs = [](double s) { return 1.0 / (s + 1.0); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinRel(std::exp(-t), TOL_STEHFEST));
    }
}

TEST_CASE("Verification func4: Talbot inverts 1/(s+1) at t=1..10",
          "[verification][func4][talbot]")
{
    nilt::Talbot algo;
    auto Fs = [](C s) -> C { return 1.0 / (s + 1.0); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinRel(std::exp(-t), TOL_TALBOT));
    }
}

TEST_CASE("Verification func4: DeHoog inverts 1/(s+1) at t=1..10",
          "[verification][func4][dehoog]")
{
    nilt::DeHoog algo;
    auto Fs = [](C s) -> C { return 1.0 / (s + 1.0); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinRel(std::exp(-t), TOL_DEHOOG));
    }
}

// func6: f(t) = t,  F(s) = 1/s^2

TEST_CASE("Verification func6: Stehfest inverts 1/s^2 at t=1..10",
          "[verification][func6][stehfest]")
{
    nilt::Stehfest algo;
    auto Fs = [](double s) { return 1.0 / (s * s); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinRel(t, TOL_STEHFEST));
    }
}

// func8: f(t) = sin(t),  F(s) = 1/(s^2+1)

TEST_CASE("Verification func8: Talbot inverts 1/(s^2+1) to sin(t) at t=1..10",
          "[verification][func8][talbot]")
{
    nilt::Talbot algo;
    auto Fs = [](C s) -> C { return 1.0 / (s * s + 1.0); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinAbs(std::sin(t), TOL_TALBOT));
    }
}

TEST_CASE("Verification func8: DeHoog inverts 1/(s^2+1) to sin(t) at t=1..10",
          "[verification][func8][dehoog]")
{
    nilt::DeHoog algo;
    auto Fs = [](C s) -> C { return 1.0 / (s * s + 1.0); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinAbs(std::sin(t), TOL_DEHOOG));
    }
}

// func10: f(t) = exp(-t)*sin(t),  F(s) = 1/((s+1)^2+1)

TEST_CASE("Verification func10: Talbot inverts damped sinusoid at t=1..10",
          "[verification][func10][talbot]")
{
    nilt::Talbot algo;
    auto Fs = [](C s) -> C { return 1.0 / ((s + 1.0) * (s + 1.0) + 1.0); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        double expected = std::exp(-t) * std::sin(t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinAbs(expected, TOL_TALBOT));
    }
}

TEST_CASE("Verification func10: DeHoog inverts damped sinusoid at t=1..10",
          "[verification][func10][dehoog]")
{
    nilt::DeHoog algo;
    auto Fs = [](C s) -> C { return 1.0 / ((s + 1.0) * (s + 1.0) + 1.0); };

    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        CAPTURE(t);
        double expected = std::exp(-t) * std::sin(t);
        REQUIRE_THAT(nilt::invert(algo, Fs, t), WithinAbs(expected, TOL_DEHOOG));
    }
}


// nilt::pi constant


TEST_CASE("nilt::pi matches std::acos(-1)", "[constants]")
{
    REQUIRE_THAT(nilt::pi, WithinRel(std::acos(-1.0), 1e-15));
}


// nilt::invert free function template deduction


TEST_CASE("nilt::invert deduces Stehfest from function pointer",
          "[api][template_deduction]")
{
    nilt::Stehfest algo;
    auto fp = +[](double s) -> double { return 1.0 / (s + 1.0); };
    double result = nilt::invert(algo, fp, 1.0);
    REQUIRE_THAT(result, WithinRel(std::exp(-1.0), 1e-4));
}
