# sage-semialgebraic

`sage-semialgebraic` is a small SageMath Python package for exact symbolic
semialgebraic sets.  Basic sets retain rational Sage polynomials, convex joins
and Cartesian products produce immutable first-order formula trees, and
coverage questions use three-valued certified results.

## Installation and quick start

Use the existing Sage Conda environment:

```bash
eval "$(conda shell.bash hook)"
conda activate sage
python -c "from sage.all import QQ; print(QQ)"
python -m pip install -e .
```

Then, from Sage's Python:

```python
from sage.all import QQ, PolynomialRing
from sage_semialgebraic import SemialgebraicSet, CoverageStatus

R = PolynomialRing(QQ, 2, names=("x", "y"))
x, y = R.gens()
square = SemialgebraicSet.from_constraints(
    ring=R, ge=[x, y], le=[x - 1, y - 1]
)

result = square.check_cube_hitting([(QQ(1)/2, QQ(1)/2)], QQ(1)/2)
assert result.status is CoverageStatus.PROVED

graphic = square.plot_projection(
    indices=(0, 1), ranges=((0, 1), (0, 1)), samples=20
)

# Ranges are optional.  Unknown/unbounded directions use a warned finite
# visualization window; this is never a boundedness or emptiness proof.
automatic_graphic = square.plot_projection()
```

For sets `A` and `B` in the same `R^n`, `A + B` is the convex join

```text
{t*a + (1-t)*b : a in A, b in B, 0 <= t <= 1},
```

not a Minkowski sum.  `A.cartesian_product(B)` preserves the coordinates of
`A` first.

## Exact and approximate operations

- Basic membership is exact, except that a caller-provided `epsilon` may be
  used only for residuals of source equality constraints.
- The rational axis-aligned box backend exactly certifies cube hitting and
  infinity-norm delta-density.  It handles strict/open faces, boundary cells,
  unbounded intervals, empty domains, and returns exact rational witnesses.
- A nonzero polynomial equality certifies that no positive-side-length ambient
  cube is contained in a basic set.  This makes cube hitting vacuous, but does
  **not** make delta-density vacuous.
- Unsupported first-order decisions return `CoverageStatus.UNKNOWN`; no grid,
  random sample, or floating-point optimizer is promoted to a proof.
- Projection plots are explicitly approximate sampled drawings.  Explicit
  finite `ranges` are honored; otherwise simple coordinate bounds are inferred
  and a warned `default_window=(-2, 2)` fills unknown or unbounded directions.
  This window is not a proof of boundedness.
- `samples`, `max_sample_points`, `max_candidates`, and `seed` accept Python
  integers and Sage `Integer` values.  A full grid is used when affordable;
  otherwise a reproducible low-discrepancy sample covers the entire window.
- Source equalities remain exact.  Plot sampling eliminates linear equalities
  over `QQ` and isolates real roots of specialized hypersurfaces in Sage's
  algebraic real field `AA`, then rechecks all constraints exactly.
- 2D/3D coordinate projections and independent custom projection bases are
  supported.  Approximation labels use supported Sage legend/text primitives,
  so `show()` and `save()` do not receive private renderer keywords.
- Nested convex joins are sampled directly with independent simplex weights.
  This preserves exact join semantics without convexifying operands and avoids
  correlated intermediate fans. Exact operand samples are mixed numerically
  only for visualization and are reported as approximate.
- `max_sample_points` is an upper bound. The normal display cloud scales like
  `4*samples^2` (minimum 1,000), capped by that argument. Sage 3D points are
  grouped linearly instead of using the slow repeated-addition path in
  `point3d(list)`.

Sage 10.5's `sage.interfaces.qepcad` module was found during development, but
the external `qepcad` executable was not installed.  The package detects this
at runtime and currently does not implement general formula translation to
QEPCAD.

## Tests

The tests are compatible with both pytest and the standard unittest runner:

```bash
eval "$(conda shell.bash hook)"
conda activate sage
python -m pytest -q
```

If pytest is not installed in the Sage environment:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

General systems containing several nonlinear equalities are not guaranteed to
have every connected component sampled.  A warning or empty drawing from the
visualizer never asserts that the exact semialgebraic set is empty.

The nested-join benchmark and measured before/after results are documented in
[`docs/PERFORMANCE.md`](docs/PERFORMANCE.md):

```bash
python benchmarks/benchmark_bicone.py --save-png
```

See `docs/MATH_SPEC.md`, `docs/ARCHITECTURE.md`, and `docs/TEST_PLAN.md` for
the precise contract and known limits.

