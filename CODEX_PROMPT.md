# Task for Codex: implement the SageMath package

You are working inside this repository on Ubuntu running under WSL. SageMath is already installed in the Conda environment named `sage` and is launched using `conda activate sage`, then `sage`. Implement a useful, tested SageMath Python package named `sage_semialgebraic`.

Read **all of** `AGENTS.md`, `docs/MATH_SPEC.md`, `docs/ARCHITECTURE.md`, `docs/TEST_PLAN.md` and `PLANS.md` before changing code. Do not stop at proposing an implementation or writing more documentation: create the Python modules, implement the algorithms and add executable pytest tests.

## Hard mathematical requirements

1. A basic set in R^n is specified by a finite conjunction of polynomial equations and inequalities over rational coefficients. Preserve exact polynomial data, support `<`, `<=`, `=`, `>=`, `>` through constructor buckets, and validate dimension/ring.
2. The operator `A + B` is precisely `{t*a+(1-t)*b : a in A, b in B, t in [0,1]}` for A,B in the same R^n. It is not Minkowski sum, and should not be replaced by convex hull of `A union B`. Store the result exactly as a quantified semialgebraic formula; do not approximate the set by a point cloud.
3. Support `A.cartesian_product(B)` in dimension n+m without variable collisions.
4. Implement 2D/3D orthogonal projection plotting using either coordinate indices or an arbitrary independent 2- or 3-vector basis. Return Sage Graphics/Graphics3d, allow explicit bounded plot ranges, and label numerical drawings as approximations. Plotting should work for representative basic sets and illustrative joins/products. Do not call interactive windows in tests.
5. Implement `check_cube_hitting(points, delta, epsilon=0, ...)`. Its meaning is: for *every* possible translated, closed, axis-aligned n-cube of side delta contained entirely in exact A, at least one accepted input point lies in that cube (including boundary). Epsilon applies to polynomial equalities when validating the input points; not to the exact definition of A or structural equality constraints.
6. General decision is first-order real arithmetic and can be intractable. For arbitrary quantified formulas use an optional certified QE/CAD backend (investigate Sage's optional `sage.interfaces.qepcad` interface). Without it, return `UNKNOWN` with a helpful explanation rather than claiming a proof based on finite grid or random sampling.
7. Regardless of QEPCAD availability, implement a **certified rational axis-aligned box fast path** for cube hitting. Suggested algorithm: reduce each admissible cube to its lower corner c; form the exact region C of allowed c. Each valid sample p hits a cube iff `c_i in [p_i-delta, p_i]` for every i. Decide whether finite union of these center boxes covers C by decomposing each coordinate at all relevant rational endpoints, testing both open intervals and singleton boundary cells. Handle open faces, unbounded boxes, empty sample sets and resource caps. Supply exact counterexample cube origins when disproved.
8. If a basic set contains a nontrivial polynomial equation, it has empty interior, hence no positive-size cube; cube hitting is `PROVED` vacuously. Separately implement `check_delta_dense`, with infinity norm and distinct semantics `for every x in A there is an accepted p with max_i |x_i-p_i| <= delta`. This is the meaningful sampling property for algebraic manifolds. Support certified elementary cases; otherwise use QE or `UNKNOWN`.
9. Return `CoverageResult` with `CoverageStatus` (`PROVED`, `DISPROVED`, `UNKNOWN`) and diagnostics: used backend, reason, and a verified witness whenever available.
10. Handle zero and constant polynomials, empty sets, strict inequalities, singleton sets, arbitrary n, open/closed boundaries, invalid delta/epsilon, dimension mismatch, input-point validation, exact vs approximate coordinates, and missing optional backends.

## Implementation priorities

- First, verify that the environment's active Python can import `sage.all`; use the existing Conda environment, do not install another copy of SageMath.
- Implement the immutable formula AST and basic API with correct existential quantification and hygienic renaming. Keep the Sage polynomial ring as the core.
- Implement simple exact `contains` fast paths, join/product and the certified box cube-hitting backend before investing heavily in general QE.
- Implement usable plotting for 2D/3D examples and a clear fallback for cases that cannot reliably be plotted.
- Add the optional QEPCAD adapter only after introspecting the actual installed Sage interface; QEPCAD may not be present. Do not introduce required dependencies or pretend it is always installed.
- Keep interfaces stable and update `README.md`, `PLANS.md` and docstrings with real implemented status.

## Test-driven acceptance criteria

Create and run tests for: join of `{0}` with `{2}` yielding segment `[0,2]` (rather than singleton `{2}`); dimension of Cartesian product; square `[0,1]^2` with delta `1/2` and point `(1/2,1/2)` returning `PROVED`; same square with only point `(0,0)` returning `DISPROVED` and a correct counterexample; nonzero equality giving vacuous hitting; unbounded boxes and finite sample points; 0/negative delta; approximate source-equality validation; custom-basis projection; WSL-headless 2D/3D plot smoke tests. Add a regression test to prevent false `PROVED` from checking only grid-aligned cubes.

Run the strongest feasible tests using the same Conda Sage interpreter, e.g. `python -m pytest -q` (or `sage -python -m pytest -q`). If something is unimplemented, either implement it or make the limitation explicit with `UNKNOWN` / a documented exception; do not silently skip core requirements or mark unverified work as complete.

Your final answer must summarize: implemented API; major exact algorithms; how to use the package from Sage; exact test commands and actual results; whether QEPCAD was detected; remaining mathematically difficult cases. Do not ask me to approve routine choices; make sensible decisions consistent with the specifications and proceed.
