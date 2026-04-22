/*
 * pybind11 bindings for the nilt (Numerical Inverse Laplace Transform) library.
 *
 * Exposes three algorithm classes (Stehfest, Talbot, DeHoog) and a free
 * function `invert(algo, F, t)` that accepts scalar or array t values.
 *
 * Stehfest: F(s) takes a real float and returns a real float.
 * Talbot/DeHoog: F(s) takes a complex and returns a complex.
 */
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>

#include "nilt.hpp"

namespace py = pybind11;

// Scalar inversion helpers
static double invert_stehfest(const nilt::Stehfest& algo,
                              py::function Fs, double t)
{
    return algo([&](double s) { return Fs(s).cast<double>(); }, t);
}

static double invert_talbot(const nilt::Talbot& algo,
                            py::function Fs, double t)
{
    return algo([&](std::complex<double> s) {
        return Fs(s).cast<std::complex<double>>();
    }, t);
}

static double invert_dehoog(const nilt::DeHoog& algo,
                            py::function Fs, double t)
{
    return algo([&](std::complex<double> s) {
        return Fs(s).cast<std::complex<double>>();
    }, t);
}

// Vectorised inversion helpers (accept numpy arrays)
static py::array_t<double> invert_stehfest_array(
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

static py::array_t<double> invert_talbot_array(
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

static py::array_t<double> invert_dehoog_array(
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

// Dispatch: picks the right scalar or array overload for any algorithm
static py::object invert_dispatch(py::object algo, py::function Fs, py::object t)
{
    // Determine algorithm type
    bool is_stehfest = py::isinstance<nilt::Stehfest>(algo);
    bool is_talbot   = py::isinstance<nilt::Talbot>(algo);
    bool is_dehoog   = py::isinstance<nilt::DeHoog>(algo);

    if (!is_stehfest && !is_talbot && !is_dehoog)
        throw std::invalid_argument("algo must be Stehfest, Talbot, or DeHoog");

    // Check if t is an ndarray
    if (py::isinstance<py::array>(t))
    {
        auto t_arr = t.cast<py::array_t<double>>();
        if (is_stehfest)
            return py::object(invert_stehfest_array(algo.cast<nilt::Stehfest>(), Fs, t_arr));
        else if (is_talbot)
            return py::object(invert_talbot_array(algo.cast<nilt::Talbot>(), Fs, t_arr));
        else
            return py::object(invert_dehoog_array(algo.cast<nilt::DeHoog>(), Fs, t_arr));
    }
    else
    {
        double t_val = t.cast<double>();
        if (is_stehfest)
            return py::cast(invert_stehfest(algo.cast<nilt::Stehfest>(), Fs, t_val));
        else if (is_talbot)
            return py::cast(invert_talbot(algo.cast<nilt::Talbot>(), Fs, t_val));
        else
            return py::cast(invert_dehoog(algo.cast<nilt::DeHoog>(), Fs, t_val));
    }
}

// Module definition
PYBIND11_MODULE(_nilt, m)
{
    m.doc() = "nilt - Numerical Inverse Laplace Transform (C++ accelerated)";

    py::class_<nilt::Stehfest>(m, "Stehfest",
        "Gaver-Stehfest algorithm for numerical inverse Laplace transform.\n"
        "F(s) must accept a real float and return a real float.")
        .def(py::init<>())
        .def_readwrite("N", &nilt::Stehfest::N,
            "Number of terms (must be even, default 18)")
        .def("__call__", &invert_stehfest,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at scalar time t")
        .def("__call__", &invert_stehfest_array,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at array of times t");

    py::class_<nilt::Talbot>(m, "Talbot",
        "Fixed Talbot algorithm for numerical inverse Laplace transform.\n"
        "F(s) must accept a complex and return a complex.")
        .def(py::init<>())
        .def_readwrite("n", &nilt::Talbot::n,
            "Number of quadrature points (default 50)")
        .def_readwrite("shift", &nilt::Talbot::shift,
            "Contour shift parameter (default 0.0)")
        .def("__call__", &invert_talbot,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at scalar time t")
        .def("__call__", &invert_talbot_array,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at array of times t");

    py::class_<nilt::DeHoog>(m, "DeHoog",
        "De Hoog et al. algorithm for numerical inverse Laplace transform.\n"
        "F(s) must accept a complex and return a complex.")
        .def(py::init<>())
        .def_readwrite("M", &nilt::DeHoog::M,
            "Order of approximation (default 40)")
        .def_readwrite("T_factor", &nilt::DeHoog::T_factor,
            "Period factor (default 4.0)")
        .def_readwrite("tol", &nilt::DeHoog::tol,
            "Tolerance (default 1e-16)")
        .def("__call__", &invert_dehoog,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at scalar time t")
        .def("__call__", &invert_dehoog_array,
            py::arg("Fs"), py::arg("t"),
            "Invert F(s) at array of times t");

    m.def("invert", &invert_dispatch,
        py::arg("algo"), py::arg("Fs"), py::arg("t"),
        "Invert the Laplace transform F(s) at time(s) t.\n\n"
        "Parameters\n"
        "----------\n"
        "algo : Stehfest, Talbot, or DeHoog\n"
        "    Algorithm instance (parameters can be tuned before calling).\n"
        "Fs : callable\n"
        "    Laplace-domain function. For Stehfest: Fs(s: float) -> float.\n"
        "    For Talbot/DeHoog: Fs(s: complex) -> complex.\n"
        "t : float or numpy.ndarray\n"
        "    Time(s) at which to evaluate the inverse. Must be positive.\n\n"
        "Returns\n"
        "-------\n"
        "float or numpy.ndarray\n"
        "    The inverse Laplace transform f(t).");

    m.attr("pi") = nilt::pi;
}
