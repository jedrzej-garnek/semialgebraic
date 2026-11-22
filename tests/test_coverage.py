from __future__ import annotations

import unittest

from sage.all import QQ, PolynomialRing

from sage_semialgebraic import CoverageStatus, SemialgebraicSet


def unit_square():
    ring = PolynomialRing(QQ, 2, names=("x", "y"))
    x, y = ring.gens()
    return SemialgebraicSet.from_constraints(
        ring=ring, ge=[x, y], le=[x - 1, y - 1])


class CoverageTests(unittest.TestCase):
    def test_square_center_hits_every_half_cube(self):
        result = unit_square().check_cube_hitting(
            [(QQ(1) / 2, QQ(1) / 2)], QQ(1) / 2)
        self.assertIs(result.status, CoverageStatus.PROVED)
        self.assertEqual(result.backend, "rational-axis-box")

    def test_square_origin_does_not_hit_every_half_cube(self):
        result = unit_square().check_cube_hitting([(0, 0)], QQ(1) / 2)
        self.assertIs(result.status, CoverageStatus.DISPROVED)
        self.assertTrue(result.metadata["witness_verified"])
        self.assertEqual(result.used_backend, "rational-axis-box")
        self.assertEqual(result.counterexample, result.witness)
        c = result.witness
        self.assertTrue(0 <= c[0] <= QQ(1) / 2)
        self.assertTrue(0 <= c[1] <= QQ(1) / 2)
        self.assertFalse(c[0] <= 0 <= c[0] + QQ(1) / 2 and
                         c[1] <= 0 <= c[1] + QQ(1) / 2)

    def test_regression_grid_aligned_origins_are_not_enough(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        interval = SemialgebraicSet.from_constraints(ring=ring, ge=[x], le=[x - 1])
        result = interval.check_cube_hitting([(0,), (1,)], QQ(1) / 2)
        self.assertIs(result.status, CoverageStatus.DISPROVED)
        c = result.witness[0]
        self.assertTrue(0 < c < QQ(1) / 2)

    def test_nonzero_equality_is_vacuously_cube_hit(self):
        ring = PolynomialRing(QQ, 2, names=("x", "y"))
        x, y = ring.gens()
        circle = SemialgebraicSet.from_constraints(ring=ring, eq=[x**2 + y**2 - 1])
        result = circle.check_cube_hitting([], QQ(1) / 100)
        self.assertIs(result.status, CoverageStatus.PROVED)
        self.assertEqual(result.backend, "empty-interior-certificate")
        self.assertTrue(result.metadata["vacuous"])

    def test_zero_equality_does_not_create_false_vacuity(self):
        ring = PolynomialRing(QQ, "x")
        whole = SemialgebraicSet.from_constraints(ring=ring, eq=[ring(0)])
        result = whole.check_cube_hitting([], 1)
        self.assertIs(result.status, CoverageStatus.DISPROVED)

    def test_unbounded_box_cannot_be_covered_by_finite_samples(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        ray = SemialgebraicSet.from_constraints(ring=ring, ge=[x])
        result = ray.check_cube_hitting([(0,), (10,)], 1)
        self.assertIs(result.status, CoverageStatus.DISPROVED)
        c = result.witness[0]
        self.assertGreaterEqual(c, 0)
        self.assertFalse(c <= 0 <= c + 1)
        self.assertFalse(c <= 10 <= c + 1)

    def test_empty_and_too_small_cube_origin_regions_are_vacuous(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        open_interval = SemialgebraicSet.from_constraints(ring=ring, gt=[x], lt=[x - 1])
        result = open_interval.check_cube_hitting([], 1)
        self.assertIs(result.status, CoverageStatus.PROVED)
        impossible = SemialgebraicSet.from_constraints(ring=ring, gt=[ring(0)])
        self.assertIs(impossible.check_cube_hitting([], 1).status,
                      CoverageStatus.PROVED)

    def test_strict_boundary_decomposition(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        open_interval = SemialgebraicSet.from_constraints(ring=ring, gt=[x], lt=[x - 1])
        result = open_interval.check_cube_hitting([(QQ(1) / 2,)], QQ(1) / 2)
        self.assertIs(result.status, CoverageStatus.PROVED)

    def test_invalid_delta_epsilon_and_point_dimension(self):
        square = unit_square()
        for delta in (0, -1):
            with self.assertRaises(ValueError):
                square.check_cube_hitting([], delta)
        with self.assertRaises(ValueError):
            square.check_cube_hitting([], 1, epsilon=-1)
        with self.assertRaises(ValueError):
            square.check_cube_hitting([(0,)], 1)

    def test_resource_cap_yields_unknown_not_guess(self):
        result = unit_square().check_cube_hitting(
            [(QQ(1) / 2, QQ(1) / 2)], QQ(1) / 2, max_cells=1)
        self.assertIs(result.status, CoverageStatus.UNKNOWN)

    def test_general_formula_without_qe_is_unknown(self):
        ring = PolynomialRing(QQ, 2, names=("x", "y"))
        x, y = ring.gens()
        disk = SemialgebraicSet.from_constraints(ring=ring, ge=[1 - x**2 - y**2])
        result = disk.check_cube_hitting([], QQ(1) / 10)
        self.assertIs(result.status, CoverageStatus.UNKNOWN)

    def test_delta_dense_box_and_singleton(self):
        square = unit_square()
        proved = square.check_delta_dense([(QQ(1) / 2, QQ(1) / 2)], QQ(1) / 2)
        self.assertIs(proved.status, CoverageStatus.PROVED)
        disproved = square.check_delta_dense([(0, 0)], QQ(1) / 2)
        self.assertIs(disproved.status, CoverageStatus.DISPROVED)

        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        singleton = SemialgebraicSet.from_constraints(ring=ring, eq=[x])
        dense = singleton.check_delta_dense([(0,)], QQ(1) / 100)
        self.assertIs(dense.status, CoverageStatus.PROVED)

    def test_product_of_boxes_uses_exact_box_backend(self):
        r = PolynomialRing(QQ, "x")
        x = r.gen()
        interval = SemialgebraicSet.from_constraints(ring=r, ge=[x], le=[x - 1])
        square = interval.cartesian_product(interval)
        result = square.check_cube_hitting([(QQ(1) / 2,) * 2], QQ(1) / 2)
        self.assertIs(result.status, CoverageStatus.PROVED)
        self.assertEqual(result.backend, "rational-axis-box")

    def test_delta_dense_manifold_is_unknown_not_vacuous(self):
        ring = PolynomialRing(QQ, 2, names=("x", "y"))
        x, y = ring.gens()
        circle = SemialgebraicSet.from_constraints(ring=ring, eq=[x**2 + y**2 - 1])
        result = circle.check_delta_dense([], 1)
        self.assertIs(result.status, CoverageStatus.UNKNOWN)

    def test_outside_samples_are_rejected(self):
        result = unit_square().check_cube_hitting([(2, 2)], QQ(1) / 2)
        self.assertIs(result.status, CoverageStatus.DISPROVED)
        self.assertEqual(result.metadata["accepted_points"], 0)
        self.assertEqual(result.metadata["rejected_points"], 1)


if __name__ == "__main__":
    unittest.main()

