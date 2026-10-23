# Package review notes

Initial review of the `reservoirflow` package on 2026-10-23. These are investigation and improvement candidates, not approved behavior changes. Confirm the intended API and numerical behavior before implementing fixes.

## Priority findings

### 1. Avoid dense matrix history in sparse simulations

**Evidence:** `reservoirflow/solutions/numerical/fdm_solution.py` initializes `self.As` as a dense array, converts sparse `A` to dense after solving (`A = A.toarray()`), then appends the flattened matrix to `self.As` on each updated step. Package search found no reads of `self.As` or `self.ds` outside this write path.

**Impact:** Retaining every matrix can consume O(nsteps × n²) memory and converting each sparse matrix to dense undermines sparse mode.

**Review direction:** Determine whether matrix history is a supported public feature. If not, remove the unused storage path; if it is, make it opt-in and preserve sparse storage where appropriate. Add a regression test that confirms sparse solves do not unexpectedly allocate dense matrix history.

### 2. Make solution clearing safe before compilation

**Evidence:** `reservoirflow/models/model.py`, `Model.clear_solutions()` deletes `self.solve` and `self.run` unconditionally. Calling this method on a new model raises `AttributeError` because those instance attributes have not yet been created. This was reproduced during review.

**Impact:** A reset/cleanup operation fails in a valid lifecycle state and cannot reliably be called repeatedly.

**Review direction:** Make clearing idempotent and ensure the model returns to its documented uncompiled state. Cover a new model, a compiled model, and repeated calls in tests.

### 3. Bound and define cache lifetimes

**Evidence:** `reservoirflow/utils/helpers.py` defines `_lru_cache` as a wrapper over `functools.lru_cache`. It is applied with `maxsize=None` to methods in `grids/regular_cartesian_grid.py`, `models/black_oil_model.py`, and `solutions/numerical/fdm_solution.py`.

**Impact:** Method cache keys include `self`, so cached entries can retain instances and results. Unbounded caches may grow with varied arguments or repeated model construction. Cached results may also become stale if inputs represented by `self` are mutated.

**Review direction:** Measure memory behavior and identify each cache's mutation dependencies. Prefer bounded/per-instance caches or explicit invalidation where needed; do not change cache policy without checking that repeated simulation calculations remain correct.

### 4. Make supported solver methods match compiler dispatch

**Evidence:** `reservoirflow/solutions/compiler.py` advertises/dispatches multiple methods whose implementation classes raise `NotImplementedError` in their constructors or solve paths. Examples include `FvmSolution` and `FemSolution`. `reservoirflow/solutions/numerical/solvers.py:get_dsolver()` is an empty `pass` stub.

**Impact:** Users can select methods that appear available but fail at runtime, and the direct-solver helper has no working behavior.

**Review direction:** Decide which methods are intended to be usable now. Either implement them or make unsupported methods unavailable with clear, early errors. Implement or remove the unused stub as appropriate, and add compiler-dispatch tests for supported and unsupported choices.

### 5. Replace bare exception handlers with specific errors

**Evidence:** Bare `except:` blocks appear at `reservoirflow/solutions/numerical/fdm_solution.py:211,500` and `reservoirflow/grids/regular_cartesian_grid.py:1055,1171`.

**Impact:** These handlers can catch unrelated programming errors and replace useful tracebacks with generic messages.

**Review direction:** Identify the expected failure at each site, catch only that exception, and chain it when raising a more helpful domain-specific error. Test both the expected invalid input and propagation of unrelated failures.

## Documentation improvements

The contribution guide (`docs/source/research_development/contribution.md`) requests NumPy-style docstrings. A repository AST scan found undocumented public methods and placeholder descriptions, especially in `BlackOilModel` and `RegularCartesianGrid`.

Prioritize docstrings for public model/grid configuration and simulation methods. Clarify:

- Input units, accepted value forms, defaults, and validation behavior.
- Return values, array shapes, indexing conventions, and whether results include boundary cells.
- Solver options, supported methods, sparse/vectorized/threading behavior, and convergence failures.
- State changes and when a compiled solution or cached result is updated.

Replace `_summary_` and `_description_` placeholders with meaningful descriptions, and correct parameter names that do not match the actual signature.

## Structure and maintainability

`reservoirflow/grids/regular_cartesian_grid.py` and `reservoirflow/models/black_oil_model.py` have accumulated broad responsibilities. For example, `BlackOilModel` includes model setup, data-frame export, and plotting delegation; the grid module contains geometry, indexing, property access, and boundary operations.

Consider extracting cohesive internal helpers/modules while preserving current public imports and behavior. Establish characterization tests before moving numerical or indexing logic, since these methods encode important shape and boundary conventions.

## Test and warning follow-up

The full test run during this review completed with **137 passed** and **54 warnings**. Warnings included:

- `reservoirflow/scalers/minmax_scaler.py:100`: `RuntimeWarning: invalid value encountered in divide` in scaler fitting.
- `reservoirflow/models/black_oil_model.py:152`: expected or defaulted initial-pressure warning in several simulation tests.

Determine whether the scaler warning represents a supported constant-range input case and define the intended transform behavior. Document or explicitly test the initial-pressure default so the warning is intentional and its effect is clear.
