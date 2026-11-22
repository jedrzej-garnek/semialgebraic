# Visualization performance

Measurements below were made on 2026-11-22 with SageMath 10.5, Python 3.11,
inside the Conda environment `sage`, with the repository on `/mnt/c` under
WSL. They are observations, not cross-machine guarantees.

Run the benchmark with:

```bash
python benchmarks/benchmark_bicone.py --save-png
```

It constructs the exact example `D = A + B + C` from the performance task.

## Confirmed baseline bottlenecks

Before optimization, warm-process exact sampling produced:

| object | points | seconds |
|---|---:|---:|
| A | 1 | 0.0084 |
| B | 407 | 0.0729 |
| C | 1 | 0.0004 |
| A+B | 50,000 | 1.4476 |
| D | 50,000 | 1.8904 |

Repeated isolated sampling of D took 3.62--4.12 seconds depending on cache and
`samples`; peak RSS for an isolated 50,000-point sample was 232,384 KiB.
`cProfile` showed exact `QQ` simplex weights, coordinate-wise Python `sum`,
deduplication, and repeated substitutions as the dominant join costs.

The larger problem occurred after sampling. Sage 10.5 implements
`point3d(list_of_points)` as repeated addition of individual `Point` graphics.
Measured construction times were 4.87 seconds for 5,000 points and 19.44
seconds for 10,000 points. The complete 50,000-point `plot_projection()` did
not finish within a controlled 30-second run.

## Optimized design

- Nested joins are flattened using associativity, without convexifying any
  individual operand.
- Singleton leaves are recognized exactly and do not consume sampling budget.
- Nonsingleton join operands receive independent low-discrepancy samples, so a
  disk or surface does not become a small collection of radial fans.
- Leaf samples remain exact and are checked against exact constraints.
- Only visualization mixtures use NumPy `float64`; `SamplingResult.exact` is
  false and a diagnostic records this boundary.
- Simplex weights and operand choices are generated in independent vectorized
  batches. Most samples fill the simplex; a smaller share covers faces and
  near-boundary regions.
- Affine equalities and remaining constraints are substituted once into the
  parameter ring instead of invoking ambient `contains()` for every candidate.
- Numerical rows are deduplicated with vectorized NumPy operations.
- 3D point clouds use one `Graphics3dGroup`, avoiding Sage's quadratic
  repeated-addition path.
- `max_sample_points` is a ceiling. The normal target is
  `min(max_sample_points, max(1000, 4*samples^2))`; `samples=50` therefore
  renders 10,000 well-distributed points instead of forcing 50,000.

The exact formula AST, exact membership, cube coverage, and delta-density code
do not use this floating-point visualization path.

## Final measurements

The complete user call with `samples=50` produced 10,000 points:

| stage | seconds |
|---|---:|
| `D.plot_projection(...)` | 0.74--0.92 |
| Tachyon PNG save (10,000 points) | 0.96--1.44 |

Peak RSS of an isolated complete call was 238,380 KiB. A deliberate
50,000-point internal sample took about 1.26 seconds; its measured stages were
1.28 seconds sampling, 0.027 seconds projection/conversion, and 1.11 seconds
for linear-time construction of the 50,000-point graphics group. Saving all
50,000 primitives through Tachyon remains expensive (33.23 seconds) and is
renderer work rather than sampling work.

Scaling with the normal display target:

| samples | target/output points | sampling seconds |
|---:|---:|---:|
| 15 | 1,000 | 0.053 |
| 30 | 3,600 | 0.224 |
| 50 | 10,000 | 0.445 |
| 100 | 40,000 | 1.242 |

For the 50,000-point geometry audit there were zero violations of the
independent bicone inequality, 19,641 positive-half points, 19,510
negative-half points, 4,998 exact middle-section operand points, and all 24
angular bins were occupied.

