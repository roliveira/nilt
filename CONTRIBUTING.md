# How to contribute

Bug reports, documentation fixes, new methods and examples, or code changes.


## Reporting issues

- Bugs: include OS, compiler or Python version, `nilt` version, and a minimal example to reproduce it.
- Feature requests: describe the use case and link to relevant papers.
- New methods or examples: see below.


## Development setup

CMake handles C++; [uv](https://docs.astral.sh/uv/) handles Python.

```bash
git clone https://github.com/roliveira/nilt.git
cd nilt

# Python (with dev deps)
uv sync --extra dev

# C++ (optional)
cmake -S . -B build -DNILT_BUILD_TESTS=ON -DNILT_BUILD_EXAMPLES=ON
cmake --build build --parallel
```


## Running tests

CI covers Ubuntu, macOS, and Windows (ARM included), Python 3.10-3.13. Everything must pass before merge.


### C++

```bash
ctest --test-dir build --output-on-failure
```


### Python

```bash
uv run pytest -v
```


## Commit messages

[Conventional Commits](https://www.conventionalcommits.org/), enforced by [python-semantic-release](https://python-semantic-release.readthedocs.io/). Every commit to `main` is checked for version bumps:

| Prefix | Effect |
|--------|--------|
| `feat` | minor bump |
| `fix`, `perf` | patch bump |
| `refactor`, `build`, `ci`, `docs`, `style`, `test`, `chore` | no bump |

```
feat: add Weeks method for numerical inversion
fix: correct tolerance in DeHoog
docs: add heat conduction example
```


## Pull requests

1. Fork and branch from `main`.
2. Follow the commit conventions.
3. Run both test suites locally.
4. Open a PR. CI runs the full matrix.
5. After merge, a release is created automatically if the commits call for one.


## Adding a new method

This is a small library, so a new method should include:

- C++ implementation with pybind11 bindings
- Tests: `tests/test_<method>.cpp` and `tests/test_<method>.py`
- Entries in the [verification suite](https://github.com/roliveira/nilt/tree/main/examples/verification)
- README updates (parameters, usage)


## Adding examples

Put new examples in the appropriate `examples/` subdirectory, document in its `README.md`, and confirm they run:

```bash
uv sync --extra plot
uv run python examples/<your_example>.py
```


## Releases

On push to `main`, CI runs tests, then `python-semantic-release` checks commits since the last tag. If a bump is warranted it updates versions, tags, and publishes a GitHub Release. You don't touch versions or tags yourself.


## GitHub Actions workflows

| Workflow | Trigger | What it does |
|----------|---------|--------------|
| `cpp-tests.yml` | PR to `main`, manual | C++ test suite |
| `python-tests.yml` | PR to `main`, manual | pytest across OS/Python matrix |
| `release.yml` | Push to `main`, manual | Tests + semantic release |
| `draft-pdf.yml` | PR to `main`, manual | JOSS paper PDF |
