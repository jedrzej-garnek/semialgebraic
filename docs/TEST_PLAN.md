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

