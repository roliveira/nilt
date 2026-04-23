[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19708531.svg)](https://doi.org/10.5281/zenodo.19708531)

# Numerical Inverse Laplace Transform Methods

A numerical inversion of Laplace transforms library implemented in C++ library with Python bindings. 
This work was partly developed in [Oliveira, R. (2021)](https://doi.org/10.25560/92253).

## Quick Start

**C++**

```cpp
#include <nilt.hpp>

// "Free" function - works with any callable 
double f = nilt::invert(nilt::Talbot{}, [](auto s) { return 1.0 / (s + 1.0); }, 1.0);

// Direct algorithm call (equivalent)
nilt::DeHoog dh;
double f = dh([](auto s) { return 1.0 / (s + 1.0); }, 2.5);

// Custom parameters
nilt::Stehfest algo;
algo.N = 12;
double f = nilt::invert(algo, my_func, 1.0);
```

**Python**

```python
import numpy as np
from nilt import Stehfest, Talbot, DeHoog, invert

# "Free" function - works with any callable
f = invert(Talbot(), lambda s: 1.0 / (s + 1.0), 1.0)

# Direct algorithm call (equivalent)
dh = DeHoog()
f = dh(lambda s: 1.0 / (s + 1.0), 2.5)

# Custom parameters
algo = Stehfest()
algo.N = 12
f = invert(algo, my_func, 1.0)

# Array of times (returns numpy array)
t = np.linspace(0.1, 10, 100)
results = invert(DeHoog(), lambda s: 1.0 / (s + 1.0), t)
```

## Methods

Three algorithms are implemented:

| C++ class        | Python class | Method       | Input           | Reference |
|------------------|--------------|--------------|-----------------|-----------|
| `nilt::Stehfest` | `nilt.Stehfest`   | Stehfest     | real `F(s)`     | [Stehfest (1970)](https://doi.org/10.1145/361953.361969)    |
| `nilt::Talbot`   | `nilt.Talbot`     | Fixed Talbot | complex `F(s)`  | [Abate & Whitt (2006)](https://doi.org/10.1287/ijoc.1050.0137) |
| `nilt::DeHoog`   | `nilt.DeHoog`     | De Hoog      | complex `F(s)`  | [De Hoog et al. (1982)](https://doi.org/10.1137/0903022)   |

All algorithms accept any callable via the free function or direct call:

|               | C++                        | Python                    |
|---------------|----------------------------|---------------------------|
| Free function | `nilt::invert(algo, F, t)` | `nilt.invert(algo, F, t)` |
| Direct call   | `algo(F, t)`               | `algo(F, t)`              |

### Parameters

Each algorithm exposes tunable parameters (identical names in C++ and Python):

| Class | Parameter | Default | Description |
|-------|-----------|---------|-------------|
| Stehfest | `N` | 18 | Number of terms (must be even) |
| Talbot   | `n` | 50 | Number of quadrature points |
| Talbot   | `shift` | 0.0 | Contour shift parameter |
| DeHoog   | `M` | 40 | Order of approximation |
| DeHoog   | `T_factor` | 4.0 | Period factor ($T = T_{\text{factor}} \cdot t$) |
| DeHoog   | `tol` | 1e-16 | Tolerance for integration limit |

## Test Functions

The [verification suite](examples/verification/) evaluates all methods against known analytical Laplace transform functions:

| # | f(t) | F(s) | Source |
|---|------|------|--------|
| 1 | $e^{-t}$ | $1/(s+1)$ | Standard |
| 2 | $1/\sqrt{\pi t}$ | $1/\sqrt{s}$ | Stehfest (1970) |
| 3 | $-\gamma - \ln t$ | $\ln(s)/s$ | Stehfest (1970) |
| 4 | $t^3/6$ | $s^{-4}$ | Stehfest (1970) |
| 5 | $\sin\sqrt{2t}$ | $\sqrt{\pi/(2s^3)}\,e^{-1/(2s)}$ | Stehfest (1970) |
| 6 | $t$ | $1/s^2$ | Abate & Whitt |
| 7 | $t\,e^{-t}$ | $1/(s+1)^2$ | Abate & Whitt |
| 8 | $\sin t$ | $1/(s^2+1)$ | Abate & Whitt |
| 9 | $\cos t$ | $s/(s^2+1)$ | Abate & Whitt |
| 10 | $e^{-t}\sin t$ | $1/((s+1)^2+1)$ | Abate & Whitt |

## Benchmark Results

See the [verification example](examples/verification/) for full the results. The table
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

```bash
cmake -B build
cmake --build build
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

Once installed, `from nilt import ...` works as expected. The `invert` function
accepts both scalar `float` and NumPy array arguments. 
Using NumPy arrays is slightly more efficient than having to evaluate several individual floats at a time. 

### Python tests (pytest)

```bash
uv run pytest                  # or simply pytest (with venv activated)
```

## Running the Verification Suite

```bash
# C++
cd examples/verification/build
./verification                 # writes CSVs to cwd
python ../plot_verification.py # reads from build/, writes PNGs there

# Python (from repo root, with .venv activated)
python examples/verification/verification.py   # writes py_*.csv to build/
```

## Examples

Several physics examples are organized by domain in `examples/`, each comparing
all three inversion methods against the known analytical solution:

| Directory       | Example              | Physics                                                   | Dimension            |
|-----------------|----------------------|-----------------------------------------------------------|----------------------|
| `verification/` | `verification`       | 10 standard test functions (Stehfest & Abate-Whitt)       | -                    |
| `transport/`    | `sphere_diffusion`   | Average concentration in a diffusing sphere               | 1D (radial)          |
| `transport/`    | `cylinder_diffusion` | Average concentration in a diffusing cylinder             | 2D (axisymmetric)    |
| `transport/`    | `advection_plume_2d` | Instantaneous release in uniform flow                     | 2D (x, y)            |
| `groundwater/`  | `theis_well`         | Drawdown from a pumping well (Theis 1935)                 | 1D (time & distance) |
| `groundwater/`  | `well_dipole`        | Pumping + injection well dipole                           | 2D (x, y)            |

Each subdirectory contains a `README.md` with the mathematical formulation and
a `plot_<example>.py` script to visualize the results. Every C++ example has a matching
Python script (`.py`) that produces identical results. Binaries are placed in a `build/`
subdirectory next to their sources; the output CSVs and PNGs are also there.


## References

- Rodolfo Oliveira. 2021. *Modelling of reactive transport in porous media using continuous time random walks*. PhD Thesis (Mar. 2021). https://doi.org/10.25560/92253
- Harald Stehfest. 1970. *Algorithm 368: Numerical inversion of Laplace transforms [D5]*. Commun. ACM 13, 1 (Jan. 1970), 47-49. https://doi.org/10.1145/361953.361969
- F.R. de Hoog, J.H. Knight, A.N. Stokes. 1982. *An improved method for numerical inversion of Laplace transforms*. SIAM J. Sci. Stat. Comput. 3, 3, 357-366. https://doi.org/10.1137/0903022
- J. Abate, W. Whitt. 2006. *A unified framework for numerically inverting Laplace transforms*. INFORMS J. Comput. 18, 4, 408-421. https://doi.org/10.1287/ijoc.1050.0137
