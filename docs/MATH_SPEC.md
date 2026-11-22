# Mathematical specification

## Basic sets

A basic set in `R^n` is a finite conjunction of relations `f rel 0`, where
`f` belongs to a Sage polynomial ring over `QQ` and `rel` is one of
`=, >=, >, <=, <`.  Constant and zero polynomials have their ordinary logical
meaning.  The ambient dimension is the number of ring generators.

An input-membership tolerance `epsilon >= 0` applies only to equalities given
to `from_constraints`: `|f(p)| <= epsilon` is accepted.  It never changes the
exact set, inequality constraints, or structural equations created by set
operations.

## Exact set operations

For `A,B subset R^n`, the overloaded addition is the convex join

```text
A + B = {t*a + (1-t)*b : a in A, b in B, 0 <= t <= 1}.
```

It is neither Minkowski addition nor, in general, `conv(A union B)`.  It is
stored as an existential formula with hygienically renamed variables.

For `A subset R^n` and `B subset R^m`, `A.cartesian_product(B)` is a subset of
`R^(n+m)` with `A` coordinates first and `B` coordinates second.

## Cube hitting

For `delta > 0`, `check_cube_hitting(P, delta)` asks whether

```text
for every c: [c,c+delta]^n subset A implies
               there exists p in P intersect A with p in [c,c+delta]^n.
```

All cubes are closed, axis parallel, and arbitrarily translated.  A finite
lattice of origins cannot certify this statement.

For an axis-aligned interval product, admissible cube origins form another
exact interval product.  A sample `p` hits precisely the origins in
`product_i [p_i-delta,p_i]`.  Coverage is decided by splitting every coordinate
at every rational endpoint and considering both open interval cells and
singleton endpoint cells.  Cartesian products of these cells have constant
membership in every sample-origin box, so checking one exact representative
per cell is a proof.

A nonzero polynomial has a zero set with empty Euclidean interior.  Therefore
a basic set containing a nonzero equality cannot contain a positive-size
ambient cube, and cube hitting is vacuously true.

## Delta-density

`check_delta_dense(P, delta, norm='linf')` has the distinct meaning

```text
for every x in A, there exists accepted p in P with
max_i |x_i-p_i| <= delta.
```

It is not vacuous on lower-dimensional algebraic sets.  The implemented exact
backend covers rational axis-aligned interval products, including singletons.

## Certification policy

`PROVED` and `DISPROVED` require an exact certificate.  An exact uncovered
rational point is returned for box-backend counterexamples.  Unsupported
general formulas return `UNKNOWN` unless a certified decision backend exists;
plots and numerical sampling are never proofs.

## Visualization sampling

Visualization uses a finite ambient window.  The window may be explicit,
inferred from independent linear coordinate bounds, or completed with a
documented default such as `[-2,2]`; none of these choices proves boundedness.
When a regular grid exceeds the candidate budget, a reproducible
low-discrepancy sample spans the entire window instead of taking a Cartesian
product prefix.

For algebraic sets, linear equalities are eliminated exactly over `QQ` and a
remaining hypersurface is specialized to univariate polynomials whose real
roots are computed in `AA`.  Returned basic-set samples satisfy all source
equalities and inequalities exactly.  Failure to find samples, especially for
several nonlinear equalities, is explicitly not an emptiness result.

For visualization, an associative tree of joins may be flattened to direct
simplex weights. This is exactly the same set operation and does not replace
any operand by its convex hull. Operand points are sampled and validated
exactly; large collections of convex combinations may then be formed in
floating point solely for plotting and must be labelled approximate.

