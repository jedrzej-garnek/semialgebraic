from __future__ import annotations

import unittest

from sage.all import QQ, PolynomialRing
from sage.plot.graphics import Graphics
from sage.plot.plot3d.base import Graphics3d

from sage_semialgebraic import SemialgebraicSet


def unit_cube():
    ring = PolynomialRing(QQ, 3, names=("x", "y", "z"))
    x, y, z = ring.gens()
    return SemialgebraicSet.from_constraints(
        ring=ring, ge=[x, y, z], le=[x - 1, y - 1, z - 1])


class PlottingTests(unittest.TestCase):
    def test_headless_coordinate_2d_plot(self):
        cube = unit_cube()
        graphic = cube.plot_projection(
            indices=(0, 2), ranges=((0, 1), (0, 1), (0, 1)), samples=3)
        self.assertIsInstance(graphic, Graphics)
        self.assertIn("approximate", graphic._extra_kwds["sage_semialgebraic_note"])

    def test_headless_coordinate_3d_plot(self):
        graphic = unit_cube().plot_projection(
            indices=(0, 1, 2), ranges=((0, 1),) * 3, samples=3)
        self.assertIsInstance(graphic, Graphics3d)

    def test_custom_basis_projection(self):
        graphic = unit_cube().plot_projection(
            basis=((1, 1, 0), (0, 0, 1)), ranges=((0, 1),) * 3, samples=3)
        self.assertIsInstance(graphic, Graphics)

    def test_dependent_basis_and_missing_ranges_are_rejected(self):
        cube = unit_cube()
        with self.assertRaises(ValueError):
            cube.plot_projection(basis=((1, 0, 0), (2, 0, 0)),
                                 ranges=((0, 1),) * 3)
        with self.assertRaises(ValueError):
            cube.plot_projection(indices=(0, 1))

    def test_join_and_product_plot_smoke(self):
        ring = PolynomialRing(QQ, "x")
        x = ring.gen()
        a = SemialgebraicSet.from_constraints(ring=ring, eq=[x])
        b = SemialgebraicSet.from_constraints(ring=ring, eq=[x - 1])
        square_as_product = (a + b).cartesian_product(a + b)
        graphic = square_as_product.plot_projection(
            indices=(0, 1), ranges=((0, 1), (0, 1)), samples=3)
        self.assertIsInstance(graphic, Graphics)


if __name__ == "__main__":
    unittest.main()

