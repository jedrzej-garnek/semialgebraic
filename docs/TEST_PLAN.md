# Test plan

The executable suite in `tests/` covers:

- every constructor relation bucket, constant/zero constraints, dimensions,
  and approximate source-equality membership;
- `{0}+{2} = [0,2]`, explicitly ruling out Minkowski behavior;
- product dimension, coordinate order, and nested variable hygiene;
- the unit-square half-cube positive and negative examples, with verified
  exact counterexample origins;
- the regression where points at `0` and `1` hit grid-aligned half-intervals
  but miss a translated interval;
- strict boundaries, empty origin regions, unbounded boxes, rejected outside
  samples, resource caps, and invalid input;
- equality-based vacuous cube hitting and non-vacuous delta-density;
- `UNKNOWN` for unsupported nonlinear coverage without a certified backend;
- coordinate/custom-basis projections and headless Sage 2D/3D graphics,
  including illustrative join/product plots.
- Sage `Integer` plotting/resource limits and rejection of nonintegral values;
- oversized-grid coverage of the whole window, reproducibility, and uniqueness;
- shifted positive-octant plotting with `samples=50` and PNG rendering;
- exact sphere sampling, multiple linear equalities, automatic windows, nested
  convex joins, 3D PNG rendering, and absence of private renderer keywords.
- independent bicone-inequality validation for flattened joins, both axial
  halves, middle-section radii, angular coverage, reproducible seeds, bounded
  output for five-way joins, and linear-time 3D point-cloud saving.

Performance measurements are separate from unit tests and live in
`benchmarks/benchmark_bicone.py`; methodology and results are recorded in
`docs/PERFORMANCE.md`.

Primary command:

```bash
python -m pytest -q
```

Dependency-free fallback inside Sage's Python:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
```

Conditional QEPCAD translation tests should be added when that backend is
implemented and the external executable is available.  They must not be
silently treated as passing when only Sage's interface module exists.

