"""Reproducible performance benchmark for nested convex joins.

Run inside the Sage Conda environment from the repository root:

    python benchmarks/benchmark_bicone.py
"""

from __future__ import annotations

import argparse
import json
import math
import resource
import tempfile
from pathlib import Path
from time import perf_counter

from sage.all import QQ, PolynomialRing

from sage_semialgebraic import SemialgebraicSet
from sage_semialgebraic.plotting import _graphics3d_point_cloud
from sage_semialgebraic.sampling import sample_set


def build_example():
    ring = PolynomialRing(QQ, names=("x", "y", "z"))
    x, y, z = ring.gens()
    top = SemialgebraicSet.from_constraints(ring=ring, eq=[x - 1, y - 1, z - 1])
    bottom = SemialgebraicSet.from_constraints(
        ring=ring, eq=[x + 1, y + 1, z + 1])
    disk = SemialgebraicSet.from_constraints(
        ring=ring, eq=[x + y + z],
        ge=[QQ(3) / 4 + x * y + y * z + z * x])
    return top, disk, bottom, top + disk + bottom


def timed_sample(set_, ranges, samples, candidates, points, seed=0):
    start = perf_counter()
    result = sample_set(set_, ranges, samples, candidates, points, seed)
    return result, perf_counter() - start


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--save-png", action="store_true")
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    top, disk, bottom, bicone = build_example()
    ranges = ((QQ(-2), QQ(2)),) * 3
    report = {"samples": [], "limits": [], "components": {}}

    for name, set_ in (("A", top), ("B", disk), ("C", bottom),
                       ("A_plus_B", top + disk), ("D", bicone)):
        result, elapsed = timed_sample(set_, ranges, 50, 50_000, 50_000)
        report["components"][name] = {
            "seconds": elapsed, "points": len(result.points)}

    for samples in (15, 30, 50, 100):
        target = min(50_000, max(1_000, 4 * samples * samples))
        result, elapsed = timed_sample(
            bicone, ranges, samples, 50_000, target)
        report["samples"].append({
            "samples": samples, "seconds": elapsed,
            "target": target, "points": len(result.points)})

    for candidates, points in ((5_000, 10_000), (10_000, 10_000),
                               (50_000, 5_000), (50_000, 10_000),
                               (50_000, 25_000), (50_000, 50_000)):
        result, elapsed = timed_sample(
            bicone, ranges, 50, candidates, points)
        report["limits"].append({
            "max_candidates": candidates, "max_sample_points": points,
            "seconds": elapsed, "points": len(result.points)})

    sampled, sample_seconds = timed_sample(
        bicone, ranges, 50, 50_000, 50_000)
    start = perf_counter()
    projected = [tuple(float(value) for value in point_)
                 for point_ in sampled.points]
    projection_seconds = perf_counter() - start
    start = perf_counter()
    staged_graphic = _graphics3d_point_cloud(
        projected, point_size=8, color="steelblue", kwds={})
    graphics_seconds = perf_counter() - start
    report["stages"] = {
        "sampling_seconds": sample_seconds,
        "projection_seconds": projection_seconds,
        "graphics3d_seconds": graphics_seconds,
    }

    violations = positive = negative = 0
    angular_bins = set()
    middle_points = 0
    for x, y, z in projected:
        coordinate_sum = x + y + z
        radius_squared = x*x + y*y + z*z - coordinate_sum**2 / 3.0
        bound = 1.5 * (1.0 - abs(coordinate_sum) / 3.0) ** 2
        violations += (abs(coordinate_sum) > 3.0 + 1e-12 or
                       radius_squared > bound + 1e-10)
        positive += coordinate_sum > 0.2
        negative += coordinate_sum < -0.2
        u = (x - y) / math.sqrt(2.0)
        v = (x + y - 2.0*z) / math.sqrt(6.0)
        radius = math.hypot(u, v)
        if radius > 0.15:
            angular_bins.add(int(
                (math.atan2(v, u) + math.pi) / (2.0*math.pi) * 24) % 24)
        middle_points += abs(coordinate_sum) < 1e-12
    report["geometry"] = {
        "violations": violations,
        "positive_half_points": positive,
        "negative_half_points": negative,
        "angular_bins_of_24": len(angular_bins),
        "exact_middle_points": middle_points,
    }

    start = perf_counter()
    graphic = bicone.plot_projection(
        indices=(0, 1, 2), ranges=ranges, samples=50)
    report["plot_projection"] = {"seconds": perf_counter() - start}
    report["max_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

    if args.save_png:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bicone.png"
            start = perf_counter()
            staged_graphic.save(str(path), viewer="tachyon")
            report["png"] = {
                "seconds": perf_counter() - start,
                "bytes": path.stat().st_size,
            }

    rendered = json.dumps(report, indent=2, sort_keys=True)
    print(rendered)
    if args.json:
        args.json.write_text(rendered + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

