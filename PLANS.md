# Implementation plan and status

- [x] Package metadata and public imports.
- [x] Immutable first-order formula AST over Sage `QQ` polynomials.
- [x] Basic constructor, relation validation, and membership with equality-only
  input tolerance.
- [x] Exact existential convex join with hygienic variables.
- [x] Exact Cartesian product with preserved coordinate order.
- [x] Certified rational axis-aligned box extraction and endpoint-cell union
  coverage, including strict faces, unbounded domains, witnesses, and caps.
- [x] Empty-interior cube-hitting certificates for nonzero equalities.
- [x] Separate infinity-norm delta-density semantics and elementary backend.
- [x] Approximate bounded 2D/3D coordinate and custom-basis Sage plotting.
- [x] Sage/Python integer validation without truncating rational inputs.
- [x] Reproducible whole-window low-discrepancy sampling when a full grid is
  over budget, with separate candidate and displayed-point limits.
- [x] Automatic finite plotting windows with coordinate-bound inference and
  explicit warnings for clipped/default directions.
- [x] Exact visualization sampling for linear equality systems and polynomial
  hypersurfaces through Sage real-algebraic root isolation.
- [x] Renderer-safe approximation labels, 2D/3D `show()` smoke checks, and PNG
  export without private `Graphics._extra_kwds` metadata.
- [x] Flatten associative join trees for visualization without convexifying
  leaves; recognize singleton operands and diversify nonsingleton samples.
- [x] Vectorize approximate join mixtures and numerical deduplication with
  NumPy while keeping exact basic samples and proof APIs separate.
- [x] Replace Sage's quadratic `point3d(list)` path with a linear
  `Graphics3dGroup` construction and scale display size predictably with
  `samples` under `max_sample_points`.
- [x] Add the bicone correctness/performance benchmark, independent geometry
  checks, measured before/after report, and saved PNG artifact.
- [x] Executable acceptance/regression tests and user documentation.
- [x] Inspect installed Sage QEPCAD interface and detect the external executable
  independently.
- [ ] Implement and validate a capture-safe AST-to-QEPCAD translator for general
  membership and coverage formulae when an executable is available.
- [ ] Add further certified backends (for example polyhedral cells) without
  weakening the `UNKNOWN` policy.

