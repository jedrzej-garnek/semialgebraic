# AGENTS.md - repository instructions for Codex

## Mission
Build a mathematical Python package for SageMath that represents semialgebraic subsets of R^n and supports the exact operation specified in `docs/MATH_SPEC.md`. This repository initially contains a scaffold, not a finished implementation. Implement working code and tests, not just further planning documents.

## Instruction precedence
- `docs/MATH_SPEC.md` is authoritative for mathematics and edge cases.
- `docs/ARCHITECTURE.md` defines the intended API and engineering boundaries.
- `docs/TEST_PLAN.md` lists acceptance criteria.
- `PLANS.md` lists implementation phases; update it as work is completed.
- If requirements conflict, preserve mathematical correctness, document the conflict, and prefer explicit `UNKNOWN` to an unproved Boolean claim.

## Environment: Ubuntu on WSL, SageMath via Conda
- Work in Linux/WSL, preferably in a project directory under the Linux home directory.
- The user has an existing Conda environment named `sage`. Do not install a second SageMath distribution or switch interpreters silently.
- Before running code, verify `conda activate sage` and `python -c "from sage.all import QQ; print(QQ)"`. In noninteractive shells, initialize Conda with `eval "$(conda shell.bash hook)"` if needed.
- Install the local editable package with `python -m pip install -e .` *after* activating the Sage environment.
- Run tests with `python -m pytest -q` or `sage -python -m pytest -q` in that same environment.
- If SageMath/pytest/QEPCAD is unavailable, report it accurately; do not claim tests passed. Never silently pip-install SageMath into a standard Python environment.
- Do not install optional native solvers or modify the user's Conda configuration without a clear need and approval.

## Immutable mathematical rules
1. For sets A,B in the same R^n, `A + B` means `{t*a + (1-t)*b : a in A, b in B, 0 <= t <= 1}`. This is a **convex join of two sets**, NOT a Minkowski sum and NOT necessarily `conv(A union B)`.
2. `A.cartesian_product(B)` lives in R^(n+m), preserving coordinates of A first, B second.
3. Mathematical set representations must be exact, symbolic and closed under quantifiers; never replace a result by a sample cloud or bounding box.
4. Arbitrary formulas need Boolean nodes and quantifiers. An AST is mandatory; a flat list of inequalities is insufficient after `+`.
5. `check_cube_hitting` means **every translated, axis-parallel, closed cube** of side `delta>0` that lies entirely in the exact set contains at least one accepted point. Checking only a finite lattice of cube origins cannot certify this claim.
6. `epsilon` permits residual error only for user-provided polynomial equalities when validating input sample points. It does not enlarge the exact set for containment of cubes, nor relax structural equations introduced by joins/products.
7. A lower-dimensional basic algebraic set defined by a nonzero polynomial equality cannot contain a positive-side-length n-dimensional cube: the hitting property is vacuously true. Provide `check_delta_dense` as a separate, non-vacuous notion.
8. Never present approximate numerical plotting or heuristic sampling as a mathematical proof.
9. Exact statements require a certified algorithm/backend. Otherwise return `UNKNOWN` (or raise a documented unsupported-operation error for APIs that require a Boolean), never a guess.

## Public interface to implement
- `SemialgebraicSet.from_constraints(ring=R, eq=(), ge=(), gt=(), le=(), lt=())`
- `ambient_dimension`, `__add__`, `cartesian_product`, `contains(point, epsilon=0)`
- `plot_projection(indices=(...), ... )` or `plot_projection(basis=(...), ...)` to 2D/3D
- `check_cube_hitting(points, delta, epsilon=0, *, method='auto', max_cells=...) -> CoverageResult`
- `check_delta_dense(points, delta, epsilon=0, *, norm='linf', method='auto', ...) -> CoverageResult`
- `CoverageStatus.PROVED`, `.DISPROVED`, `.UNKNOWN`, and explainable result metadata.
See `docs/ARCHITECTURE.md` for detailed signatures and feasible backend choices.

## Engineering and testing rules
- Keep public code importable with `from sage_semialgebraic import SemialgebraicSet, CoverageStatus`.
- Use Sage polynomial rings, `QQ`, algebraic tools and graphics; avoid SymPy as a competing core representation.
- Prefer pure Python `.py` modules importing `sage.all` over `.sage` preparser-dependent source.
- Be careful with distinct polynomial rings, ambient dimensions, variable renaming and quantifier capture.
- Use exact rational arithmetic for certified box algorithms; do not convert `QQ` to `float` before proof steps.
- Require explicit finite drawing windows when needed; do not assume unbounded sets fit a picture.
- Keep plotting usable headlessly in WSL; do not call GUI `show()` in unit tests.
- Validate `delta > 0`, `epsilon >= 0`, point dimension, independent basis vectors, and equal dimensions for `+`.
- Add fast paths for basic-set membership, products, singleton joins, zero-interior equalities, and exact rational axis-aligned boxes.
- QEPCAD through Sage is an optional adapter, not a mandatory dependency. Detect it at runtime and add conditional tests.
- Maintain type hints, API docstrings with math definitions, robust errors and deterministic tests. Avoid unrelated dependencies.
- After any implementation step, run the full feasible test suite, describe what was actually tested and keep `README.md` and `PLANS.md` synchronized.

## Completion criteria
- A clean editable installation inside Conda `sage`.
- Public operations return mathematically correct symbolic sets.
- Plotting works on at least representative 2D/3D examples and clearly identifies heuristic output.
- Cube hitting is certified for axis-aligned rational boxes, handles empty/interiorless cases, and reports `UNKNOWN` for unsupported general formulas unless a certified QE backend is available.
- No fabricated all-cubes guarantees from grids, sampling or floating-point optimizers.
- Tests cover the corner cases listed in `docs/TEST_PLAN.md`.
- Final Codex report: files changed; runnable examples; tests and outcomes; unsupported cases; exact next steps.
