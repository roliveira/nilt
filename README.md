[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19708530.svg)](https://doi.org/10.5281/zenodo.19708530) [![status](https://joss.theoj.org/papers/01af96de50bd5e2aa96b6658243f9b9c/status.svg)](https://joss.theoj.org/papers/01af96de50bd5e2aa96b6658243f9b9c)

# NILT: Numerical Inverse Laplace Transform Methods

A C++ header-only library (with Python bindings) for numerically inverting Laplace transforms[^1]. Three algorithms are provided: Gaver-Stehfest, fixed Talbot, and De Hoog et al. All share the same callable interface in both languages.

[^1]: This work was partly developed in [Oliveira, R. (2021)](https://doi.org/10.25560/92253).


## Statement of need

Many problems in physics and engineering are easier to solve in the Laplace domain than in the time domain. Groundwater drawdown, heat conduction in semi-infinite solids, diffusion from spheres and cylinders, viscoelastic creep are great examples that have closed-form Laplace-domain solutions that are difficult or impossible to invert analytically.

Existing tools are scattered:

- MATLAB's `ilaplace` implements an inverse Laplace transform but it has no access to individual methods or parameters within it, and does not offer an open-source license.
- Python's `mpmath.invertlaplace` provides all three families of methods (and Cohen method as well) and is written in pure Python with arbitrary-precision arithmetic, so a Python-first implementation is far slower when you need to invert at thousands of points.
- The [`ilt`](https://github.com/nocliper/ilt) package wraps a single algorithm and it provides an implementation that is too tightly integrated to the application (transient spectroscopy). 
- No other C++ library packages multiple algorithms behind a common interface.

NILT provides Stehfest, Talbot, and De Hoog in a dependency-free C++ header that compiles with any C++14 toolchain. The Python bindings expose the same compiled code for scripting and prototyping. 


## Quick Start

**C++**

```cpp
#include <nilt.hpp>

// "Free" function - works with any callable 
double f = nilt::invert(nilt::Talbot{}, [](auto s) { return 1.0 / (s + 1.0); }, 1.0);

// Direct algorithm call (equivalent)
nilt::DeHoog dh;
double f = dh([](auto s) { return 1.0 / (s + 1.0); }, 2.5);

// Vector of times (algorithm is constructed once, reused for all t)
std::vector<double> t = {0.1, 0.5, 1.0, 2.0, 5.0};
auto results = nilt::invert(nilt::Talbot{}, [](auto s) { return 1.0 / (s + 1.0); }, t);

// Custom parameters (see Parameters section for full list)
nilt::Stehfest algo;
algo.N = 12;
double f = nilt::invert(algo, my_func, 1.0);
```

**Python**

```python
import nilt

# Scalar evaluation (Stehfest is the default method)
f = nilt.invert(lambda s: 1.0 / (s + 1.0), 1.0)

# Array of times (returns numpy array)
f = nilt.invert(lambda s: 1.0 / (s + 1.0), [0.1, 0.5, 1.0, 2.0, 5.0])

# Pick a different method, pass parameters as keyword arguments
f = nilt.invert(lambda s: 1.0 / (s + 1.0), 1.0, method="Talbot", N=64)

# Class instances are callable directly (useful when reusing an algorithm)
algo = nilt.DeHoog()
algo.M = 60
f = algo(lambda s: 1.0 / (s + 1.0), 2.5)
```


## Methods

Three algorithms are implemented:

| C++ class        | Python class | Method       | Input           | Reference |
|------------------|--------------|--------------|-----------------|-----------|
| `nilt::Stehfest` | `nilt.Stehfest`   | Stehfest     | real `F(s)`     | [Stehfest (1970)](https://doi.org/10.1145/361953.361969)    |
| `nilt::Talbot`   | `nilt.Talbot`     | Fixed Talbot | complex `F(s)`  | [Abate & Whitt (2006)](https://doi.org/10.1287/ijoc.1050.0137) |
| `nilt::DeHoog`   | `nilt.DeHoog`     | De Hoog      | complex `F(s)`  | [De Hoog et al. (1982)](https://doi.org/10.1137/0903022)   |

All algorithms accept any callable via the free function or direct call:

|               | C++                        | Python                                  |
|---------------|----------------------------|------------------------------------------|
| Free function | `nilt::invert(algo, F, t)` | `nilt.invert(F, t, method=..., **kwargs)` |
| Direct call   | `algo(F, t)`               | `algo(F, t)`                             |

The Python free function selects the algorithm by name (case-insensitive) and
forwards keyword arguments as parameters. The C++ free function takes an
algorithm instance directly.


### Parameters

Each algorithm exposes tunable parameters (identical names in C++ and Python):

| Class | Parameter | Default | Description |
|-------|-----------|---------|-------------|
| Stehfest | `N` | 18 | Number of terms (must be even) |
| Talbot   | `N` | 50 | Number of quadrature points |
| Talbot   | `SHIFT` | 0.0 | Contour shift parameter |
| DeHoog   | `M` | 40 | Order of approximation |
| DeHoog   | `T_FACTOR` | 4.0 | Period factor ($T = T_{\text{FACTOR}} \cdot t$) |
| DeHoog   | `TOL` | 1e-16 | Tolerance for integration limit |


## Test Functions

The [verification suite](verification/) evaluates all methods against known analytical Laplace transform functions:

| # | f(t) | F(s) | Source |
|---|------|------|--------|
| 1 | $1/\sqrt{\pi t}$ | $1/\sqrt{s}$ | Stehfest (1970) |
| 2 | $-\gamma - \ln t$ | $\ln(s)/s$ | Stehfest (1970) |
| 3 | $t^3/6$ | $s^{-4}$ | Stehfest (1970) |
| 4 | $e^{-t}$ | $1/(s+1)$ | Standard |
| 5 | $\sin\sqrt{2t}$ | $\sqrt{\pi/(2s^3)}\,e^{-1/(2s)}$ | Stehfest (1970) |
| 6 | $t$ | $1/s^2$ | Abate & Whitt |
| 7 | $t\,e^{-t}$ | $1/(s+1)^2$ | Abate & Whitt |
| 8 | $\sin t$ | $1/(s^2+1)$ | Abate & Whitt |
| 9 | $\cos t$ | $s/(s^2+1)$ | Abate & Whitt |
| 10 | $e^{-t}\sin t$ | $1/((s+1)^2+1)$ | Abate & Whitt |


## Benchmark Results

See the [verification directory](verification/) for the full results. The table
below shows a test function from Stehfest (1970) ($f(t) = 1/\sqrt{\pi t}$) as an example:

| t | f(t) | Stehfest | err | Talbot | err | De Hoog | err |
|---|------|----------|-----|--------|-----|---------|-----|
| 1 | 5.6419e-01 | 5.6419e-01 | 2.17e-06 | 5.6419e-01 | 4.63e-12 | 5.6419e-01 | 1.73e-13 |
| 2 | 3.9894e-01 | 3.9894e-01 | 4.92e-06 | 3.9894e-01 | 4.82e-12 | 3.9894e-01 | 2.70e-14 |
| 3 | 3.2574e-01 | 3.2573e-01 | 6.34e-06 | 3.2574e-01 | 2.74e-12 | 3.2574e-01 | 2.11e-14 |
| 4 | 2.8209e-01 | 2.8210e-01 | 2.17e-06 | 2.8209e-01 | 4.63e-12 | 2.8209e-01 | 1.73e-13 |
| 5 | 2.5231e-01 | 2.5231e-01 | 4.24e-06 | 2.5231e-01 | 4.87e-12 | 2.5231e-01 | 5.06e-14 |
| 6 | 2.3033e-01 | 2.3033e-01 | 8.70e-07 | 2.3033e-01 | 2.54e-12 | 2.3033e-01 | 7.58e-14 |
| 7 | 2.1324e-01 | 2.1324e-01 | 2.81e-06 | 2.1324e-01 | 5.25e-12 | 2.1324e-01 | 4.14e-14 |
| 8 | 1.9947e-01 | 1.9947e-01 | 4.92e-06 | 1.9947e-01 | 4.82e-12 | 1.9947e-01 | 2.70e-14 |
| 9 | 1.8806e-01 | 1.8806e-01 | 6.24e-06 | 1.8806e-01 | 4.61e-12 | 1.8806e-01 | 3.26e-14 |
| 10 | 1.7841e-01 | 1.7841e-01 | 5.70e-06 | 1.7841e-01 | 4.84e-12 | 1.7841e-01 | 6.02e-14 |


## Building

### C++ library

The library is built and installed from `CMakeLists.txt` using CMake (+3.19). If you're just installing the library, make sure to turn off the examples using `-DNILT_BUILD_EXAMPLES=OFF`.

Install the headers and CMake config files to a chosen prefix:

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=/path/to/install
cmake --build build
cmake --install build
```

Then consume from another CMake project:

```cmake
find_package(nilt REQUIRED)
target_link_libraries(my_target PRIVATE nilt::nilt)
```

### C++ tests (Catch2, fetched automatically)

```bash
cmake -B build -DNILT_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```


### Python bindings (pybind11)

The Python package is built and installed automatically from `pyproject.toml`
using [scikit-build-core](https://scikit-build-core.readthedocs.io/), which
drives the CMake build behind the scenes.

With [uv](https://docs.astral.sh/uv/) (recommended):

```bash
uv sync --extra dev  # creates venv, builds C++ extension, installs everything
```

Or with pip:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Once installed, `import nilt` works as expected. The `invert` function
accepts scalars, lists, tuples, and NumPy arrays as time arguments.
Passing an array is more efficient than calling `invert` in a loop.


### Python tests (pytest)

```bash
uv run pytest                  # or simply pytest (with venv activated)
```


## Running the Verification Suite

```bash
# C++ (from repo root)
cmake -B build -DNILT_BUILD_VERIFICATION=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build
cd build
./verification
./benchmark_timing

# Python (from build/ directory)
python ../verification/verification.py
python ../verification/benchmark_timing.py
```


## Examples

Several physics examples are organized by domain in `examples/`, each comparing
all three inversion methods against the known analytical solution:

| Directory               | Example              | Physics                                                   | Dimension            |
|-------------------------|----------------------|-----------------------------------------------------------|----------------------|
| `examples/transport/`   | `sphere_diffusion`   | Average concentration in a diffusing sphere               | 1D (radial)          |
| `examples/transport/`   | `cylinder_diffusion` | Average concentration in a diffusing cylinder             | 2D (axisymmetric)    |
| `examples/transport/`   | `advection_plume_2d` | Instantaneous release in uniform flow                     | 2D (x, y)            |
| `examples/groundwater/` | `theis_well`         | Drawdown from a pumping well (Theis 1935)                 | 1D (time & distance) |
| `examples/groundwater/` | `well_dipole`        | Pumping + injection well dipole                           | 2D (x, y)            |

Each subdirectory contains a `README.md` with the mathematical formulation and
a `plot_<example>.py` script to visualize the results. Every C++ example has a
matching Python script that produces identical output. All scripts write to the
current working directory, so run them from `build/`:

```bash
cd build
./groundwater_theis_well                           # C++ -> cpp_*.csv
python ../examples/groundwater/theis_well.py       # Python -> py_*.csv
python ../examples/groundwater/plot_theis_well.py  # reads *_.csv, writes *.png
```


## Contributing

Contributions for bugs, features, other methods and examples are all welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for the development setup, commit conventions, pull request guidelines and etc.


## References

- Rodolfo Oliveira. 2021. *Modelling of reactive transport in porous media using continuous time random walks*. PhD Thesis (Mar. 2021). https://doi.org/10.25560/92253
- Harald Stehfest. 1970. *Algorithm 368: Numerical inversion of Laplace transforms [D5]*. Commun. ACM 13, 1 (Jan. 1970), 47-49. https://doi.org/10.1145/361953.361969
- F.R. de Hoog, J.H. Knight, A.N. Stokes. 1982. *An improved method for numerical inversion of Laplace transforms*. SIAM J. Sci. Stat. Comput. 3, 3, 357-366. https://doi.org/10.1137/0903022
- J. Abate, W. Whitt. 2006. *A unified framework for numerically inverting Laplace transforms*. INFORMS J. Comput. 18, 4, 408-421. https://doi.org/10.1287/ijoc.1050.0137
