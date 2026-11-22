from __future__ import annotations

import tempfile
import unittest
import warnings
from pathlib import Path

from sage.all import QQ, ZZ, PolynomialRing
from sage.plot.graphics import Graphics
from sage.plot.plot3d.base import Graphics3d

from sage_semialgebraic import Exists, SemialgebraicSet
from sage_semialgebraic.plotting import PlottingWindowWarning
from sage_semialgebraic.sampling import candidate_points, sample_set


def unit_cube():
    ring = PolynomialRing(QQ, 3, names=("x", "y", "z"))
    x, y, z = ring.gens()
    return SemialgebraicSet.from_constraints(
        ring=ring, ge=[x, y, z], le=[x - 1, y - 1, z - 1])


class PlottingTests(unittest.TestCase):
    def test_shifted_positive_octant_uses_whole_sampling_window(self):
        ring = PolynomialRing(QQ, names=("x", "y", "z"))
        x, y, z = ring.gens()
        octant = SemialgebraicSet.from_constraints(ring=ring, ge=[x - 1, y - 1, z - 1])
        self.assertTrue(octant.contains((1, 1, 1)))
        self.assertTrue(octant.contains((2, 2, 2)))
        self.assertFalse(octant.contains((0, 0, 0)))
        graphic = octant.plot_projection(
            ranges=[(-2, 2)] * 3, samples=50, max_sample_points=50_000)
        self.assertIsInstance(graphic, Graphics)
        self.assertGreater(len(graphic), 0)
        self.assertNotIn("sage_semialgebraic_note", graphic._extra_kwds)
        graphic.matplotlib()
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "octant.png"
            graphic.save(str(filename))
            self.assertGreater(filename.stat().st_size, 0)

    def test_sage_integer_limits_and_nonintegral_rejection(self):
        cube = unit_cube()
        graphic = cube.plot_projection(
            ranges=((-1, 2),) * 3, samples=ZZ(12),
            max_sample_points=ZZ(2_000), seed=ZZ(7))
        self.assertGreater(len(graphic), 0)
        for invalid in (2.5, QQ(3) / 2):
            with self.assertRaises(ValueError):
                cube.plot_projection(ranges=((0, 1),) * 3, samples=invalid)
        with self.assertRaises(ValueError):
            cube.plot_projection(ranges=((0, 1),) * 3,
                                 max_sample_points=QQ(3) / 2)

    def test_oversized_grid_is_reproducible_and_not_a_prefix(self):
        ranges = ((QQ(-2), QQ(2)),) * 3
        first = tuple(candidate_points(ranges, 50, 100, 19))
        second = tuple(candidate_points(ranges, 50, 100, 19))
        self.assertEqual(first, second)
        self.assertEqual(len(first), len(set(first)))
        self.assertTrue(any(all(value >= 1 for value in point) for point in first))
        self.assertTrue(any(point[0] > 0 and point[1] < 0 for point in first))

    def test_sphere_exact_sampling_and_2d_3d_plots(self):
        ring = PolynomialRing(QQ, names=("x", "y", "z"))
        x, y, z = ring.gens()
        sphere = SemialgebraicSet.from_constraints(ring=ring, eq=[x**2 + y**2 + z**2 - 1])
        self.assertTrue(sphere.contains((1, 0, 0)))
        self.assertTrue(sphere.contains((0, 1, 0)))
        self.assertFalse(sphere.contains((0, 0, 0)))
        ranges = ((QQ(-3) / 2, QQ(3) / 2),) * 3
        sampled = sample_set(sphere, ranges, 16, 6_000, 6_000, 0)
        self.assertGreater(len(sampled.points), 0)
        self.assertTrue(all(sphere.contains(point) for point in sampled.points[:50]))
        graphic2d = sphere.plot_projection(
            indices=(0, 1), ranges=ranges, samples=16, max_sample_points=6_000)
        graphic3d = sphere.plot_projection(
            indices=(0, 1, 2), ranges=ranges, samples=16,
            max_sample_points=6_000)
        self.assertGreater(len(graphic2d), 0)
        self.assertIsInstance(graphic3d, Graphics3d)
        with tempfile.TemporaryDirectory() as directory:
            filename = Path(directory) / "sphere.png"
            graphic3d.save(str(filename), viewer="tachyon")
            self.assertGreater(filename.stat().st_size, 0)

    def test_algebraic_curve_surface_intersection_and_cylinder(self):
        ring = PolynomialRing(QQ, names=("x", "y", "z"))
        x, y, z = ring.gens()
        ranges = ((QQ(-3) / 2, QQ(3) / 2),) * 3
        great_circle = SemialgebraicSet.from_constraints(
            ring=ring, eq=[x**2 + y**2 + z**2 - 1, z])
        circle_points = sample_set(great_circle, ranges, 12, 4_000, 4_000)
        self.assertGreater(len(circle_points.points), 0)
        self.assertTrue(all(great_circle.contains(point)
                            for point in circle_points.points[:50]))

        cylinder = SemialgebraicSet.from_constraints(
            ring=ring, eq=[x**2 + y**2 - 1], ge=[z], le=[z - 2])
        cylinder_ranges = ((-1.5, 1.5), (-1.5, 1.5), (0, 2))
        cylinder_graphic = cylinder.plot_projection(
            indices=(0, 1, 2), ranges=cylinder_ranges,
            samples=12, max_sample_points=4_000)
        self.assertIsInstance(cylinder_graphic, Graphics3d)
        lower, upper = cylinder_graphic.bounding_box()
        self.assertLess(lower[2], upper[2])

    def test_multiple_linear_equalities_sample_singleton(self):
        ring = PolynomialRing(QQ, names=("x", "y", "z"))
        x, y, z = ring.gens()
        singleton = SemialgebraicSet.from_constraints(ring=ring, eq=[x, y, z])
        self.assertTrue(singleton.contains((0, 0, 0)))
        self.assertFalse(singleton.contains((1, 0, 0)))
        sampled = sample_set(singleton, ((-2, 2),) * 3, 10, 100, 100)
        self.assertEqual(sampled.points, ((0, 0, 0),))
        graphic = singleton.plot_projection(ranges=((-2, 2),) * 3)
        self.assertGreater(len(graphic), 0)

        cone_point = SemialgebraicSet.from_constraints(
            ring=ring, eq=[x + y + z], ge=[x, y, z])
        cone_samples = sample_set(cone_point, ((-2, 2),) * 3, 10, 100, 100)
        self.assertIn((0, 0, 0), cone_samples.points)

    def test_univariate_algebraic_equality_sampling(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        roots = SemialgebraicSet.from_constraints(ring=ring, eq=[x**2 - 2])
        sampled = sample_set(roots, ((-2, 2),), 10, 100, 100)
        self.assertEqual(len(sampled.points), 2)
        self.assertTrue(all(roots.contains(point) for point in sampled.points))

    def test_automatic_ranges_and_default_projection(self):
        ring = PolynomialRing(QQ, names=("x", "y", "z"))
        x, y, z = ring.gens()
        octant = SemialgebraicSet.from_constraints(ring=ring, ge=[x - 1, y - 1, z - 1])
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            graphic = octant.plot_projection(samples=8, max_sample_points=1_000)
        self.assertGreater(len(graphic), 0)
        self.assertTrue(any(issubclass(item.category, PlottingWindowWarning)
                            for item in caught))

    def test_headless_coordinate_and_custom_basis_plots(self):
        cube = unit_cube()
        graphic2d = cube.plot_projection(
            indices=(0, 2), ranges=((0, 1),) * 3, samples=3)
        graphic3d = cube.plot_projection(
            indices=(0, 1, 2), ranges=((0, 1),) * 3, samples=3)
        custom = cube.plot_projection(
            basis=((1, 1, 0), (0, 0, 1)), ranges=((0, 1),) * 3, samples=3)
        self.assertIsInstance(graphic2d, Graphics)
        self.assertIsInstance(graphic3d, Graphics3d)
        self.assertIsInstance(custom, Graphics)
        with self.assertRaises(ValueError):
            cube.plot_projection(basis=((1, 0, 0), (2, 0, 0)),
                                 ranges=((0, 1),) * 3)
        with self.assertRaises(ValueError):
            cube.plot_projection(indices=(0, QQ(3) / 2), ranges=((0, 1),) * 3)

    def test_nested_convex_joins_of_inequality_and_equality_sets(self):
        ring = PolynomialRing(QQ, 2, names=("x", "y"))
        x, y = ring.gens()
        a = SemialgebraicSet.from_constraints(ring=ring, eq=[x, y])
        b = SemialgebraicSet.from_constraints(ring=ring, eq=[x - 1, y])
        c = SemialgebraicSet.from_constraints(
            ring=ring, ge=[x], le=[x - 1], eq=[y - 1])
        joined = a + b + c
        self.assertIsInstance(joined.formula, Exists)
        graphic = joined.plot_projection(
            ranges=((-1, 2), (-1, 2)), samples=10,
            max_sample_points=5_000, seed=3)
        self.assertGreater(len(graphic), 0)


if __name__ == "__main__":
    unittest.main()
