/*
 * Snapshot tests for numerical inverse Laplace transform methods.
 *
 * These tests lock in the exact numerical output at full double precision for
 * all 10 standard verification functions x 3 methods.  They complement the
 * tolerance-based tests: tolerance tests verify analytical correctness,
 * snapshots detect any unintended change in numerical output.
 *
 * To update approved files after intentional algorithm changes, run the tests
 * and copy the .received.txt files over the .approved.txt files (or use an
 * ApprovalTests reporter).
 */
#include <catch2/catch_test_macros.hpp>

#include "ApprovalTests/Approvals.h"
#include "ApprovalTests/namers/NamerFactory.h"

#include "nilt.hpp"
#include "testfunctions.hpp"

#include <cmath>
#include <complex>
#include <iomanip>
#include <sstream>
#include <string>

using C = std::complex<double>;

namespace {

template <typename Algo, typename FuncS>
std::string serialize(Algo& algo, FuncS Fs)
{
    std::ostringstream oss;
    oss << "t,ftn\n";
    oss << std::scientific << std::setprecision(4);
    for (int ti = 1; ti <= 10; ++ti) {
        double t = static_cast<double>(ti);
        double ftn = nilt::invert(algo, Fs, t);
        oss << std::fixed << std::setprecision(1) << t << ","
            << std::scientific << std::setprecision(4) << ftn << "\n";
    }
    return oss.str();
}

} // namespace


// Stehfest (all 10 functions, real F(s))

TEST_CASE("Snapshot Stehfest func1", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs1<double>));
}

TEST_CASE("Snapshot Stehfest func2", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs2<double>));
}

TEST_CASE("Snapshot Stehfest func3", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs3<double>));
}

TEST_CASE("Snapshot Stehfest func4", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs4<double>));
}

TEST_CASE("Snapshot Stehfest func5", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs5<double>));
}

TEST_CASE("Snapshot Stehfest func6", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs6<double>));
}

TEST_CASE("Snapshot Stehfest func7", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs7<double>));
}

TEST_CASE("Snapshot Stehfest func8", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs8<double>));
}

TEST_CASE("Snapshot Stehfest func9", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs9<double>));
}

TEST_CASE("Snapshot Stehfest func10", "[snapshot][stehfest]")
{
    nilt::Stehfest algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs10<double>));
}


// Talbot (complex F(s), skip func5 - branch cut)

TEST_CASE("Snapshot Talbot func1", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs1<C>));
}

TEST_CASE("Snapshot Talbot func2", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs2<C>));
}

TEST_CASE("Snapshot Talbot func3", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs3<C>));
}

TEST_CASE("Snapshot Talbot func4", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs4<C>));
}

TEST_CASE("Snapshot Talbot func6", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs6<C>));
}

TEST_CASE("Snapshot Talbot func7", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs7<C>));
}

TEST_CASE("Snapshot Talbot func8", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs8<C>));
}

TEST_CASE("Snapshot Talbot func9", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs9<C>));
}

TEST_CASE("Snapshot Talbot func10", "[snapshot][talbot]")
{
    nilt::Talbot algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs10<C>));
}


// DeHoog (complex F(s), skip func5 - branch cut)

TEST_CASE("Snapshot DeHoog func1", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs1<C>));
}

TEST_CASE("Snapshot DeHoog func2", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs2<C>));
}

TEST_CASE("Snapshot DeHoog func3", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs3<C>));
}

TEST_CASE("Snapshot DeHoog func4", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs4<C>));
}

TEST_CASE("Snapshot DeHoog func6", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs6<C>));
}

TEST_CASE("Snapshot DeHoog func7", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs7<C>));
}

TEST_CASE("Snapshot DeHoog func8", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs8<C>));
}

TEST_CASE("Snapshot DeHoog func9", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs9<C>));
}

TEST_CASE("Snapshot DeHoog func10", "[snapshot][dehoog]")
{
    nilt::DeHoog algo;
    ApprovalTests::Approvals::verify(serialize(algo, Fs10<C>));
}
