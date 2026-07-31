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

#include "nilt.hpp"
#include "util.hpp"


namespace py = pybind11;

// --- Stehfest helpers (real-valued F(s)) ------------------------------------

static double stehfest_scalar(const nilt::Stehfest& algo,
                              py::function Fs, double t)
{
    return algo([&](double s) { return Fs(s).cast<double>(); }, t);
}

static py::array_t<double> stehfest_array(
    const nilt::Stehfest& algo, py::function Fs, py::array_t<double> t_arr)
{
    auto t = t_arr.unchecked<1>();
    py::array_t<double> out(t.shape(0));
    auto o = out.mutable_unchecked<1>();
    auto cpp_Fs = [&](double s) { return Fs(s).cast<double>(); };
    for (py::ssize_t i = 0; i < t.shape(0); ++i)
        o(i) = algo(cpp_Fs, t(i));
    return out;
}

// --- Talbot helpers (complex-valued F(s)) -----------------------------------

static double talbot_scalar(const nilt::Talbot& algo,
                            py::function Fs, double t)
{
    return algo([&](std::complex<double> s) {
        return Fs(s).cast<std::complex<double>>();
    }, t);
}

static py::array_t<double> talbot_array(
    const nilt::Talbot& algo, py::function Fs, py::array_t<double> t_arr)
{
    auto t = t_arr.unchecked<1>();
    py::array_t<double> out(t.shape(0));
    auto o = out.mutable_unchecked<1>();
    auto cpp_Fs = [&](std::complex<double> s) {
        return Fs(s).cast<std::complex<double>>();
    };
    for (py::ssize_t i = 0; i < t.shape(0); ++i)
        o(i) = algo(cpp_Fs, t(i));
    return out;
}

// --- DeHoog helpers (complex-valued F(s)) -----------------------------------

static double dehoog_scalar(const nilt::DeHoog& algo,
                            py::function Fs, double t)
{
    return algo([&](std::complex<double> s) {
        return Fs(s).cast<std::complex<double>>();
    }, t);
}

static py::array_t<double> dehoog_array(
    const nilt::DeHoog& algo, py::function Fs, py::array_t<double> t_arr)
{
    auto t = t_arr.unchecked<1>();
    py::array_t<double> out(t.shape(0));
    auto o = out.mutable_unchecked<1>();
    auto cpp_Fs = [&](std::complex<double> s) {
        return Fs(s).cast<std::complex<double>>();
    };
    for (py::ssize_t i = 0; i < t.shape(0); ++i)
        o(i) = algo(cpp_Fs, t(i));
    return out;
}

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
