/*
 * pybind11 bindings for the nilt (Numerical Inverse Laplace Transform) library.
 *
 * Exposes three algorithm classes (Stehfest, Talbot, DeHoog) with scalar
 * and array __call__ overloads.  The high-level invert() API lives in
 * Python (nilt/__init__.py) and delegates to these classes.
 *
 * Stehfest: F(s) takes a real float and returns a real float.
 * Talbot/DeHoog: F(s) takes a complex and returns a complex.
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <cstring>

#include "nilt.hpp"
#include "util.hpp"


namespace py = pybind11;

// --- Batched forwarders ------------------------------------------------------
// Each algo exposes two `eval_batched` overloads: per-t (single s-array) and
// per-array (all s-values across all t in one shot).  Both are forwarded via
// a lambda that (1) copies the s-array into a numpy array, (2) calls Python
// Fs once, (3) copies the result back.  Result: one Python round-trip per
// scalar call and one round-trip per whole-array call.

template<typename S>
struct NumpyBatched
{
    py::function Fs;
    void operator()(const S* s, S* out, int n) const
    {
        py::array_t<S> s_arr(n);
        std::memcpy(s_arr.mutable_data(), s, n * sizeof(S));
        auto f_arr = Fs(s_arr).template cast<py::array_t<S>>();
        std::memcpy(out, f_arr.data(), n * sizeof(S));
    }
};

template<typename Algo, typename S>
static double scalar_call(const Algo& algo, py::function Fs, double t)
{
    return algo.eval_batched(NumpyBatched<S>{Fs}, t);
}

template<typename Algo, typename S>
static py::array_t<double> array_call(
    const Algo& algo, py::function Fs, py::array_t<double> t_arr)
{
    py::array_t<double> out(t_arr.shape(0));
    algo.eval_batched(
        NumpyBatched<S>{Fs},
        t_arr.data(), out.mutable_data(), static_cast<int>(t_arr.shape(0)));
    return out;
}

// --- Stehfest ---------------------------------------------------------------

static double stehfest_scalar(const nilt::Stehfest& a, py::function Fs, double t)
    { return scalar_call<nilt::Stehfest, double>(a, Fs, t); }

static py::array_t<double> stehfest_array(
    const nilt::Stehfest& a, py::function Fs, py::array_t<double> t)
    { return array_call<nilt::Stehfest, double>(a, Fs, t); }

// --- Talbot -----------------------------------------------------------------

static double talbot_scalar(const nilt::Talbot& a, py::function Fs, double t)
    { return scalar_call<nilt::Talbot, std::complex<double>>(a, Fs, t); }

static py::array_t<double> talbot_array(
    const nilt::Talbot& a, py::function Fs, py::array_t<double> t)
    { return array_call<nilt::Talbot, std::complex<double>>(a, Fs, t); }

// --- DeHoog -----------------------------------------------------------------

static double dehoog_scalar(const nilt::DeHoog& a, py::function Fs, double t)
    { return scalar_call<nilt::DeHoog, std::complex<double>>(a, Fs, t); }

static py::array_t<double> dehoog_array(
    const nilt::DeHoog& a, py::function Fs, py::array_t<double> t)
    { return array_call<nilt::DeHoog, std::complex<double>>(a, Fs, t); }

// --- Module definition ------------------------------------------------------

PYBIND11_MODULE(_nilt, m)
{
    m.doc() = "nilt - Numerical Inverse Laplace Transform (C++ accelerated)";

    py::class_<nilt::Stehfest>(m, "Stehfest",
        "Gaver-Stehfest algorithm (real-valued F(s)).")
        .def(py::init<>())
        .def_readwrite("N", &nilt::Stehfest::N,
            "Number of terms (must be even, 2-20, default 18)")
        .def("__call__", &stehfest_scalar,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at scalar time t")
        .def("__call__", &stehfest_array,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at array of times t");

    py::class_<nilt::Talbot>(m, "Talbot",
        "Fixed Talbot algorithm (complex-valued F(s)).")
        .def(py::init<>())
        .def_readwrite("N", &nilt::Talbot::N,
            "Number of quadrature points (default 50, table-accelerated 8-64)")
        .def_readwrite("SHIFT", &nilt::Talbot::SHIFT,
            "Real-axis contour shift (default 0.0)")
        .def("__call__", &talbot_scalar,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at scalar time t")
        .def("__call__", &talbot_array,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at array of times t");

    py::class_<nilt::DeHoog>(m, "DeHoog",
        "De Hoog et al. algorithm (complex-valued F(s)).")
        .def(py::init<>())
        .def_readwrite("M", &nilt::DeHoog::M,
            "Order of approximation (default 40)")
        .def_readwrite("T_FACTOR", &nilt::DeHoog::T_FACTOR,
            "Period factor T = T_FACTOR * t (default 4.0)")
        .def_readwrite("TOL", &nilt::DeHoog::TOL,
            "Contour damping tolerance (default 1e-16)")
        .def("__call__", &dehoog_scalar,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at scalar time t")
        .def("__call__", &dehoog_array,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at array of times t");

    m.attr("pi") = nilt::util::PI;
}
