from __future__ import annotations

import math
import unittest

from sage.all import QQ, PolynomialRing

from sage_semialgebraic import Exists, SemialgebraicSet


class CoreTests(unittest.TestCase):
    def test_all_relation_buckets_and_basic_membership(self):
        ring = PolynomialRing(QQ, 2, names=("x", "y"))
        x, y = ring.gens()
        set_ = SemialgebraicSet.from_constraints(
            ring=ring, eq=[x - y], ge=[x], gt=[x + 1], le=[x - 2], lt=[x - 3])
        self.assertTrue(set_.contains((1, 1)))
        self.assertFalse(set_.contains((-1, -1)))
        self.assertFalse(set_.contains((1, 0)))

    def test_join_of_singletons_is_segment_not_minkowski_sum(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        zero = SemialgebraicSet.from_constraints(ring=ring, eq=[x])
        two = SemialgebraicSet.from_constraints(ring=ring, eq=[x - 2])
        segment = zero + two
        self.assertIsInstance(segment.formula, Exists)
        self.assertTrue(segment.contains((0,)))
        self.assertTrue(segment.contains((1,)))
        self.assertTrue(segment.contains((2,)))
        self.assertFalse(segment.contains((-1,)))
        self.assertFalse(segment.contains((3,)))

    def test_join_dimension_mismatch(self):
        r1 = PolynomialRing(QQ, "x")
        r2 = PolynomialRing(QQ, 2, names=("x", "y"))
        a = SemialgebraicSet.from_constraints(ring=r1)
        b = SemialgebraicSet.from_constraints(ring=r2)
        with self.assertRaises(ValueError):
            _ = a + b

    def test_cartesian_product_dimension_and_coordinate_order(self):
        r1 = PolynomialRing(QQ, "u")
        u = r1.gen()
        r2 = PolynomialRing(QQ, 2, names=("u", "v"))
        x, y = r2.gens()
        a = SemialgebraicSet.from_constraints(ring=r1, eq=[u - 1])
        b = SemialgebraicSet.from_constraints(ring=r2, eq=[x - 2, y - 3])
        product = a.cartesian_product(b)
        self.assertEqual(product.ambient_dimension, 3)
        self.assertTrue(product.contains((1, 2, 3)))
        self.assertFalse(product.contains((2, 1, 3)))

    def test_nested_operations_have_hygienic_variables(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        a = SemialgebraicSet.from_constraints(ring=ring, eq=[x])
        b = SemialgebraicSet.from_constraints(ring=ring, eq=[x - 1])
        nested = (a + b).cartesian_product(a + b)
        self.assertEqual(nested.ambient_dimension, 2)
        self.assertTrue(nested.contains((QQ(1) / 2, QQ(3) / 4)))

    def test_approximate_source_equality_validation(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        algebraic = SemialgebraicSet.from_constraints(ring=ring, eq=[x**2 - 2])
        point = (math.sqrt(2.0),)
        self.assertFalse(algebraic.contains(point, epsilon=0))
        self.assertTrue(algebraic.contains(point, epsilon=1e-12))
        with self.assertRaises(ValueError):
            algebraic.contains(point, epsilon=-1)

    def test_constant_and_zero_polynomials(self):
        ring = PolynomialRing(QQ, "x")
        impossible = SemialgebraicSet.from_constraints(ring=ring, gt=[ring(0)])
        whole = SemialgebraicSet.from_constraints(ring=ring, eq=[ring(0)], ge=[ring(1)])
        self.assertFalse(impossible.contains((0,)))
        self.assertTrue(whole.contains((123,)))

    def test_point_dimension_is_validated(self):
        ring = PolynomialRing(QQ, 2, names=("x", "y"))
        set_ = SemialgebraicSet.from_constraints(ring=ring)
        with self.assertRaises(ValueError):
            set_.contains((0,))


if __name__ == "__main__":
    unittest.main()

