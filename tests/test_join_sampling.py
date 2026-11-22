from __future__ import annotations

import math
import tempfile
import unittest
import warnings
from pathlib import Path

from sage.all import QQ, PolynomialRing
from sage.plot.plot3d.base import Graphics3d

from sage_semialgebraic import Exists, SemialgebraicSet
from sage_semialgebraic.sampling import sample_set


def bicone_example():
    ring = PolynomialRing(QQ, names=("x", "y", "z"))
    x, y, z = ring.gens()
    top = SemialgebraicSet.from_constraints(ring=ring, eq=[x - 1, y - 1, z - 1])
    bottom = SemialgebraicSet.from_constraints(
        ring=ring, eq=[x + 1, y + 1, z + 1])
    disk = SemialgebraicSet.from_constraints(
        ring=ring, eq=[x + y + z],
        ge=[QQ(3) / 4 + x * y + y * z + z * x])
    return ring, top, disk, bottom, top + disk + bottom


class JoinSamplingTests(unittest.TestCase):
    def test_bicone_geometry_distribution_and_exact_symbolics(self):
        _, top, disk, bottom, bicone = bicone_example()
        self.assertIsInstance(bicone.formula, Exists)
        self.assertTrue(top.contains((1, 1, 1)))
        self.assertTrue(bottom.contains((-1, -1, -1)))
        self.assertTrue(disk.contains((0, 0, 0)))

        ranges = ((QQ(-2), QQ(2)),) * 3
        result = sample_set(bicone, ranges, 30, 8_000, 8_000, seed=17)
        self.assertEqual(len(result.points), 8_000)
        self.assertFalse(result.exact)

        positive = negative = 0
        angular_bins = set()
        radii = []
        middle_radii = []
        for x, y, z in result.points:
            coordinate_sum = x + y + z
            radius_squared = x * x + y * y + z * z - coordinate_sum**2 / 3.0
            bound = 1.5 * (1.0 - abs(coordinate_sum) / 3.0) ** 2
            self.assertLessEqual(abs(coordinate_sum), 3.0 + 1e-12)
            self.assertLessEqual(radius_squared, bound + 1e-10)
            positive += coordinate_sum > 0.2
            negative += coordinate_sum < -0.2
            u = (x - y) / math.sqrt(2.0)
            v = (x + y - 2.0 * z) / math.sqrt(6.0)
            radius = math.hypot(u, v)
            radii.append(radius)
            if radius > 0.15:
                angle = math.atan2(v, u) + math.pi
                angular_bins.add(int(angle / (2.0 * math.pi) * 16) % 16)
            if abs(coordinate_sum) < 1e-12:
                middle_radii.append(radius)

        self.assertGreater(positive, 1_000)
        self.assertGreater(negative, 1_000)
        self.assertGreaterEqual(len(angular_bins), 14)
        self.assertLess(min(radii), 0.02)
        self.assertGreater(max(radii), 1.0)
        self.assertGreater(len(middle_radii), 30)
        self.assertLess(min(middle_radii), 0.1)
        self.assertGreater(max(middle_radii), 1.0)

    def test_join_sampling_is_reproducible(self):
        *_, bicone = bicone_example()
        ranges = ((QQ(-2), QQ(2)),) * 3
        first = sample_set(bicone, ranges, 20, 2_000, 2_000, seed=123)
        second = sample_set(bicone, ranges, 20, 2_000, 2_000, seed=123)
        different = sample_set(bicone, ranges, 20, 2_000, 2_000, seed=124)
        self.assertEqual(first.points, second.points)
        self.assertNotEqual(first.points, different.points)

    def test_flattened_join_has_bounded_output_and_3d_save(self):
        ring, top, disk, bottom, bicone = bicone_example()
        x, y, z = ring.gens()
        extra_left = SemialgebraicSet.from_constraints(
            ring=ring, eq=[x - 1, y + 1, z])
        extra_right = SemialgebraicSet.from_constraints(
            ring=ring, eq=[x + 1, y - 1, z])
        five_way = top + disk + bottom + extra_left + extra_right
        ranges = ((QQ(-2), QQ(2)),) * 3
        result = sample_set(five_way, ranges, 20, 3_000, 3_000, seed=4)
        self.assertEqual(len(result.points), 3_000)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            graphic = bicone.plot_projection(
                indices=(0, 1, 2), ranges=ranges, samples=20,
                max_candidates=2_000, max_sample_points=2_000, seed=9)
        self.assertIsInstance(graphic, Graphics3d)
        lower, upper = graphic.bounding_box()
        self.assertLess(lower[0], upper[0])
        self.assertLess(lower[1], upper[1])
        self.assertLess(lower[2], upper[2])
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "bicone.png"
            graphic.save(str(filename), viewer="tachyon")
            self.assertGreater(filename.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()

