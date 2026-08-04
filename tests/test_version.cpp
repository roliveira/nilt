#include <catch2/catch_test_macros.hpp>

#include "nilt.hpp"

#include <cstring>
#include <string>

// Compile-time gate: if this ever fires, either the generated version.hpp is
// stale or `nilt.hpp` failed to include it.  Loose lower bound so this test
// keeps passing across 3.x -> 4.x without needing an edit each release.
static_assert(NILT_VERSION_MAJOR >= 3,
    "NILT_VERSION_MAJOR unexpectedly low; check the generated nilt/version.hpp");

TEST_CASE("Version macros are exposed", "[version]")
{
    // Non-negative components.  MINOR/PATCH can legitimately be 0.
    REQUIRE(NILT_VERSION_MAJOR >= 0);
    REQUIRE(NILT_VERSION_MINOR >= 0);
    REQUIRE(NILT_VERSION_PATCH >= 0);

    // NILT_VERSION_STRING must be a non-empty literal that starts with a digit.
    REQUIRE(std::strlen(NILT_VERSION_STRING) > 0);
    REQUIRE(NILT_VERSION_STRING[0] >= '0');
    REQUIRE(NILT_VERSION_STRING[0] <= '9');
}

TEST_CASE("Version string matches MAJOR.MINOR.PATCH", "[version]")
{
    // Rebuild the canonical string from the numeric macros and compare.
    std::string expected =
        std::to_string(NILT_VERSION_MAJOR) + "." +
        std::to_string(NILT_VERSION_MINOR) + "." +
        std::to_string(NILT_VERSION_PATCH);
    REQUIRE(std::string(NILT_VERSION_STRING) == expected);
}
