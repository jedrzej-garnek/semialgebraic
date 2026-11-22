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
set_.plot_projection(indices=(...), ranges=(...), ...)
set_.plot_projection(basis=(...), ranges=(...), ...)
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
- `plotting.py`: bounded numerical sampling and orthogonal coordinate/custom
  basis projection into Sage 2D/3D graphics.
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

Plotting is intentionally separate from certification.  It requires a finite
range for every ambient coordinate, returns Sage graphics without calling
`show()`, and labels output as an approximate sampled projection.  A custom
basis contains two or three independent ambient vectors; Gram-Schmidt produces
an orthonormal coordinate system for the span.

