#include <cmath>
#include <complex>
#include <functional>
#include <string>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include "nilt.hpp"
#include "testfunctions.hpp"
#include "conftest.hpp"

using C = std::complex<double>;

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;


TEST_CASE("Test real verification functions: Stehfest", "[testfunctions][stehfest]")
{
    int index = 1;
    for (const auto& test_case : cases_real) {
        DYNAMIC_SECTION("Evaluating Stehfest: f" << index++) {

            std::string name_str = test_case.name;
            if (name_str == "f4" || name_str == "f7" || name_str == "f8" || name_str == "f9" || name_str == "f10") {
                SKIP("Stehfest diverges for oscillatory functions"); 
            }

            for (int ti = 1; ti <= 10; ++ti) {
                double t = static_cast<double>(ti);
                CAPTURE(t);
                
                nilt::Stehfest stehfest;
                double s_val = nilt::invert(stehfest, test_case.Fs, t);
                double expected = test_case.ft(t);

                REQUIRE_THAT(s_val, WithinRel(expected, TOL["STEHFEST_REL"].get<double>()));
            }
        }
    }
}


TEST_CASE("Test complex verification functions: Talbot and DeHoog", "[testfunctions][talbot,dehoog]")
{
    int index = 1;
    for (const auto& test_case : cases_complex) {
        DYNAMIC_SECTION("Evaluating Talbot/DeHoog: f" << index++) {
            for (int ti = 1; ti <= 10; ++ti) {
                double t = static_cast<double>(ti);
                CAPTURE(t);
                
                nilt::Talbot   talbot;
                nilt::DeHoog   dehoog;

                double t_val = nilt::invert(talbot, test_case.Fs, t);
                double d_val = nilt::invert(dehoog, test_case.Fs, t);
                double expected = test_case.ft(t);

                REQUIRE_THAT(t_val, WithinRel(expected, TOL["TALBOT_REL"].get<double>()));
                REQUIRE_THAT(d_val, WithinRel(expected, TOL["DEHOOG_REL"].get<double>()));
            }
        }
    }
}
