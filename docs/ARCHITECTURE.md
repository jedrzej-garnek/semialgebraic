# Architecture

## Public API

```python
SemialgebraicSet.from_constraints(*, ring, eq=(), ge=(), gt=(), le=(), lt=())
set_.ambient_dimension
set_.formula
set_.ring
set_.contains(point, epsilon=0)
set_a + set_b
set_a.cartesian_product(set_b)
set_.check_cube_hitting(points, delta, epsilon=0,
                        method='auto', max_cells=200_000)
set_.check_delta_dense(points, delta, epsilon=0, norm='linf',
                       method='auto', max_cells=200_000)
set_.plot_projection(indices=(...), ranges=None, samples=15,
                     max_sample_points=50_000, max_candidates=None,
                     seed=0, default_window=(-2, 2), ...)
set_.plot_projection(basis=(...), ...)
```

Top-level imports include `SemialgebraicSet`, `CoverageResult`,
`CoverageStatus`, and the formula node types.

## Modules

- `formula.py`: frozen Boolean/atomic/existential AST nodes, exact polynomial
  remapping, and formatting.
- `set.py`: basic constructor, hygienic joins/products, membership fast paths,
  and public dispatch.
- `coverage.py`: rational interval recognition, accepted-point filtering,
  empty-interior certificates, endpoint-cell decomposition, witnesses, and
  QEPCAD availability detection.
- `sampling.py`: full-grid/low-discrepancy candidate generation, exact linear
  elimination, exact algebraic hypersurface specialization, recursive
  join/product sampling, vectorized numerical join mixtures, and automatic
  finite-window inference.
- `plotting.py`: orthogonal coordinate/custom-basis projection into Sage 2D/3D
  graphics and safe user-facing sampling diagnostics.
- `result.py`: immutable three-valued results and diagnostics.

Every symbolic result uses Sage polynomials over `QQ`.  A result ring contains
the free ambient generators followed by fresh internal generators.  Quantifier
nodes refer to fresh generator names; remapping recursively copies polynomial
monomials, preventing capture across nested joins and products.

## Supported decisions

Basic conjunction membership and products are implemented directly.  Join
membership has an exact singleton-to-singleton segment fast path.  Other
quantified membership raises a documented `NotImplementedError` because the
Boolean API has no `UNKNOWN` return type.

Coverage supports exact rational axis-aligned boxes, their Cartesian products,
one-dimensional joins of recognized singletons, empty sets recognized by that
backend, and equality-based empty-interior certificates.  General nonlinear
or quantified cases return `CoverageStatus.UNKNOWN`.

Sage's QEPCAD Python interface is optional.  Runtime availability requires both
that module and the external `qepcad` executable.  General AST-to-QEPCAD
translation is a future backend and is not advertised as implemented.

## Plotting contract

Plotting is intentionally separate from certification.  Explicit finite ranges
are honored; otherwise simple coordinate bounds are inferred and unbounded or
unknown directions use a configurable finite default window with a warning.
This never proves boundedness or emptiness.  Oversized grids use a reproducible
Halton sample over the whole window rather than a lexicographic prefix.

Linear equalities are eliminated over `QQ`.  A remaining polynomial
hypersurface is sampled by specializing all but one variable and isolating
exact real roots in Sage's `AA`; every returned basic-set point is then checked
with exact `contains(..., epsilon=0)`.  General systems of several nonlinear
equalities are sampled component-wise where possible and emit diagnostics;
component coverage is not guaranteed.

The method returns Sage graphics without calling `show()` itself and labels
output through supported legend/text primitives.  It never stores private
metadata in `Graphics._extra_kwds`.  A custom basis contains two or three
independent ambient vectors; Gram-Schmidt produces orthonormal coordinates for
the same span.

Nested joins are flattened because the binary join is associative, but no leaf
is replaced by its convex hull. Operand samples are exact; their visualization
mixtures are explicitly marked approximate and use NumPy `float64`. The
renderer receives one `Graphics3dGroup`, avoiding Sage 10.5's quadratic
`point3d(list)` addition path. `max_sample_points` is a ceiling; the normal
target grows as `4*samples^2` between a floor of 1,000 and that ceiling.

See `docs/PERFORMANCE.md` for the executable benchmark and measured costs.

