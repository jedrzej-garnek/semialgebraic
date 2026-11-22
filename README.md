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
- Projection plots are explicitly approximate sampled drawings.  Finite
  ambient `ranges` are required, and 2D/3D coordinate or independent custom
  bases are accepted.

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

See `docs/MATH_SPEC.md`, `docs/ARCHITECTURE.md`, and `docs/TEST_PLAN.md` for
the precise contract and known limits.

