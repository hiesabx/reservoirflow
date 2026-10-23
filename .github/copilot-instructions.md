# Copilot instructions

## Project overview

ReservoirFlow is a Python reservoir-simulation library (Python 3.11+). The installable package is `reservoirflow/`; its public entry point is `import reservoirflow as rf`, which exposes the main subpackages.

A simulation is assembled from a grid (geometry and rock properties), fluid properties, wells, and a model. These domain objects inherit shared unit, dtype, and verbosity behavior from `Base`. Models such as `BlackOilModel` combine those inputs and expose the simulation state. Calling `model.compile(stype, method)` uses `solutions/compiler.py` to select a numerical, analytical, or neurical solution implementation. The compiled solution owns `solve()`/`run()` and simulation results; the model forwards those methods and can keep multiple compiled solutions. When changing solver behavior, trace both the model setup and the selected solution implementation.

The main domain packages are `grids`, `fluids`, `wells`, and `models`; solution implementations are grouped under `solutions/numerical`, `solutions/analytical`, and `solutions/neurical`. `backends` contains array-backend adapters, while `plots`, `scalers`, and `utils` provide supporting functionality. Documentation is built from `docs/source` and includes API references and tutorial material.

## Commands

The project uses `uv` and declares pytest in the `dev` dependency group:

```console
uv sync --group dev
uv run pytest
uv run pytest tests/test_example_7_1.py
uv run pytest tests/test_example_7_1.py::test_pressures
uv build
```

Ruff is configured in `pyproject.toml` for Python 3.11; when Ruff is available, run:

```console
ruff check .
```

To build the Sphinx documentation after syncing its dependencies:

```console
uv sync --group docs
uv run sphinx-build -b html docs/source docs/build/html
```

## Repository conventions

- Follow the contribution guide at `docs/source/research_development/contribution.md`: Python style follows PEP 8/20, and docstrings use NumPy style. Use built-in collection type annotations as in PEP 585.
- Keep the package’s unit system consistent. `Base` supports `field`, `metric`, and `lab` units and carries unit-specific labels and conversion factors; related simulation components commonly accept `unit`, `dtype`, and `verbose`.
- Tests live in `tests/` and use pytest. Simulation tests commonly construct a small grid/fluid/model, compile a solver, and compare results with NumPy assertions. Sparse/dense and vectorized/symbolic behavior is tested through pytest parametrization; preserve those dimensions when changing solver paths.
- User-facing examples and docs import the library as `rf` (`import reservoirflow as rf`); keep public API examples consistent with that convention.
