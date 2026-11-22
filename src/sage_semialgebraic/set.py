"""Core exact semialgebraic-set representation."""

from __future__ import annotations

from itertools import count
from typing import Iterable, Sequence

from sage.all import QQ, PolynomialRing, vector

from .formula import (And, Atom, Exists, Formula, Relation, Truth,
                      format_formula, polynomial_variable_names, remap_formula)


_BUCKETS = ("eq", "ge", "gt", "le", "lt")
_RELATIONS = {
    "eq": Relation.EQ,
    "ge": Relation.GE,
    "gt": Relation.GT,
    "le": Relation.LE,
    "lt": Relation.LT,
}


def _fresh_ring(names: Sequence[str]):
    if not names:
        return PolynomialRing(QQ, "_unused")
    return PolynomialRing(QQ, len(names), names=tuple(names))


def _validate_epsilon(epsilon: object) -> None:
    try:
        if epsilon < 0:
            raise ValueError("epsilon must be nonnegative")
    except TypeError as exc:
        raise TypeError("epsilon must be an ordered real number") from exc


class SemialgebraicSet:
    """An exact symbolic semialgebraic subset of a real affine space.

    Use :meth:`from_constraints` to create a basic set.  Formulae produced by
    joins and products retain their existential variables in an immutable AST.
    """

    __slots__ = (
        "_ambient_dimension", "_ambient_names", "_ring", "_formula",
        "_kind", "_constraints", "_operands",
    )

    def __init__(self, *, ambient_dimension: int, ambient_names: Sequence[str],
                 ring: object, formula: Formula, kind: str,
                 constraints: dict[str, tuple[object, ...]] | None = None,
                 operands: tuple["SemialgebraicSet", ...] = ()) -> None:
        self._ambient_dimension = int(ambient_dimension)
        self._ambient_names = tuple(ambient_names)
        self._ring = ring
        self._formula = formula
        self._kind = kind
        self._constraints = constraints
        self._operands = operands

    @classmethod
    def from_constraints(cls, *, ring: object, eq: Iterable[object] = (),
                         ge: Iterable[object] = (), gt: Iterable[object] = (),
                         le: Iterable[object] = (), lt: Iterable[object] = ()) \
            -> "SemialgebraicSet":
        """Create a basic set from rational-polynomial constraints.

        Every bucket compares its polynomials with zero.  Thus ``ge=[x]``
        means ``x >= 0`` and ``le=[x-1]`` means ``x-1 <= 0``.
        """
        if not hasattr(ring, "gens") or not hasattr(ring, "base_ring"):
            raise TypeError("ring must be a Sage polynomial ring")
        if ring.base_ring() is not QQ and ring.base_ring() != QQ:
            raise TypeError("the polynomial ring must have rational coefficients")
        generators = tuple(ring.gens())
        if not generators:
            raise ValueError("the polynomial ring must have at least one generator")
        supplied = {"eq": eq, "ge": ge, "gt": gt, "le": le, "lt": lt}
        constraints: dict[str, tuple[object, ...]] = {}
        atoms: list[Formula] = []
        for bucket in _BUCKETS:
            try:
                polynomials = tuple(ring(p) for p in supplied[bucket])
            except (TypeError, ValueError) as exc:
                raise TypeError(f"a {bucket} constraint is not in the supplied ring") from exc
            constraints[bucket] = polynomials
            atoms.extend(Atom(p, _RELATIONS[bucket]) for p in polynomials)
        return cls(
            ambient_dimension=len(generators),
            ambient_names=tuple(str(g) for g in generators),
            ring=ring,
            formula=And(atoms) if atoms else Truth(True),
            kind="basic",
            constraints=constraints,
        )

    @property
    def ambient_dimension(self) -> int:
        return self._ambient_dimension

    @property
    def formula(self) -> Formula:
        """The immutable exact first-order formula defining this set."""
        return self._formula

    @property
    def ring(self) -> object:
        """Polynomial ring containing all free and quantified variables."""
        return self._ring

    def __repr__(self) -> str:
        return (f"SemialgebraicSet(R^{self.ambient_dimension}, "
                f"{format_formula(self.formula)})")

    def _embedding_plan(self, prefix: str, free_names: Sequence[str]):
        all_names = polynomial_variable_names(self.formula)
        internal = [n for n in all_names if n not in self._ambient_names]
        mapped_internal = [f"{prefix}q{i}" for i in range(len(internal))]
        name_map = dict(zip(self._ambient_names, free_names))
        name_map.update(zip(internal, mapped_internal))
        return name_map, mapped_internal

    def __add__(self, other: object) -> "SemialgebraicSet":
        """Return the exact convex join of two sets in the same ``R^n``.

        The result is ``{t*a + (1-t)*b | a in self, b in other, 0<=t<=1}``.
        This operation is deliberately not Minkowski addition.
        """
        if not isinstance(other, SemialgebraicSet):
            return NotImplemented
        if self.ambient_dimension != other.ambient_dimension:
            raise ValueError("convex join requires equal ambient dimensions")
        n = self.ambient_dimension
        x_names = [f"x{i}" for i in range(n)]
        a_names = [f"a{i}" for i in range(n)]
        b_names = [f"b{i}" for i in range(n)]
        amap, a_internal = self._embedding_plan("left_", a_names)
        bmap, b_internal = other._embedding_plan("right_", b_names)
        tau_name = "tau"
        names = x_names + a_names + b_names + a_internal + b_internal + [tau_name]
        ring = _fresh_ring(names)
        by_name = {str(g): g for g in ring.gens()}
        left = remap_formula(self.formula, ring,
                             {k: by_name[v] for k, v in amap.items()})
        right = remap_formula(other.formula, ring,
                              {k: by_name[v] for k, v in bmap.items()})
        tau = by_name[tau_name]
        structural = [
            Atom(by_name[x_names[i]] -
                 (tau * by_name[a_names[i]] +
                  (1 - tau) * by_name[b_names[i]]), Relation.EQ)
            for i in range(n)
        ]
        body = And([left, right, Atom(tau, Relation.GE),
                    Atom(1 - tau, Relation.GE), *structural])
        quantified = a_names + b_names + [tau_name]
        formula = Exists(quantified, body)
        return SemialgebraicSet(
            ambient_dimension=n,
            ambient_names=x_names,
            ring=ring,
            formula=formula,
            kind="join",
            operands=(self, other),
        )

    def cartesian_product(self, other: "SemialgebraicSet") -> "SemialgebraicSet":
        """Return ``self x other`` with self's coordinates first."""
        if not isinstance(other, SemialgebraicSet):
            raise TypeError("cartesian_product expects a SemialgebraicSet")
        n, m = self.ambient_dimension, other.ambient_dimension
        x_names = [f"x{i}" for i in range(n + m)]
        amap, a_internal = self._embedding_plan("left_", x_names[:n])
        bmap, b_internal = other._embedding_plan("right_", x_names[n:])
        names = x_names + a_internal + b_internal
        ring = _fresh_ring(names)
        by_name = {str(g): g for g in ring.gens()}
        left = remap_formula(self.formula, ring,
                             {k: by_name[v] for k, v in amap.items()})
        right = remap_formula(other.formula, ring,
                              {k: by_name[v] for k, v in bmap.items()})
        return SemialgebraicSet(
            ambient_dimension=n + m,
            ambient_names=x_names,
            ring=ring,
            formula=And([left, right]),
            kind="product",
            operands=(self, other),
        )

    def _validate_point(self, point: Sequence[object]) -> tuple[object, ...]:
        try:
            result = tuple(point)
        except TypeError as exc:
            raise TypeError("point must be an iterable") from exc
        if len(result) != self.ambient_dimension:
            raise ValueError(
                f"point has dimension {len(result)}; expected {self.ambient_dimension}")
        return result

    def contains(self, point: Sequence[object], epsilon: object = 0) -> bool:
        """Test membership when an exact implemented fast path is available.

        ``epsilon`` affects only equality constraints supplied to basic-set
        constructors.  General quantified membership raises
        :class:`NotImplementedError` rather than guessing.
        """
        _validate_epsilon(epsilon)
        p = self._validate_point(point)
        if self._kind == "basic":
            substitutions = dict(zip(self._ring.gens(), p))
            for bucket in _BUCKETS:
                for polynomial in self._constraints[bucket]:  # type: ignore[index]
                    value = polynomial.subs(substitutions)
                    if bucket == "eq":
                        if epsilon == 0:
                            ok = value == 0
                        else:
                            try:
                                ok = abs(value) <= epsilon
                            except TypeError:
                                ok = abs(float(value)) <= float(epsilon)
                    elif bucket == "ge":
                        ok = value >= 0
                    elif bucket == "gt":
                        ok = value > 0
                    elif bucket == "le":
                        ok = value <= 0
                    else:
                        ok = value < 0
                    if not bool(ok):
                        return False
            return True
        if self._kind == "product":
            left, right = self._operands
            n = left.ambient_dimension
            return left.contains(p[:n], epsilon=epsilon) and right.contains(
                p[n:], epsilon=epsilon)
        if self._kind == "join":
            from .coverage import singleton_point
            left, right = self._operands
            a, b = singleton_point(left), singleton_point(right)
            if a is not None and b is not None:
                difference = vector(QQ, b) - vector(QQ, a)
                target = vector(QQ, p) - vector(QQ, a)
                if difference.is_zero():
                    return target.is_zero()
                parameter = None
                for d, q in zip(difference, target):
                    if d != 0:
                        parameter = q / d
                        break
                assert parameter is not None
                return (0 <= parameter <= 1 and
                        target == parameter * difference)
        raise NotImplementedError(
            "membership for this quantified formula needs a certified QE backend")

    def check_cube_hitting(self, points: Iterable[Sequence[object]], delta: object,
                           epsilon: object = 0, *, method: str = "auto",
                           max_cells: int = 200_000):
        """Certify hitting of every contained closed axis-aligned delta-cube."""
        from .coverage import check_cube_hitting
        return check_cube_hitting(self, points, delta, epsilon,
                                  method=method, max_cells=max_cells)

    def check_delta_dense(self, points: Iterable[Sequence[object]], delta: object,
                          epsilon: object = 0, *, norm: str = "linf",
                          method: str = "auto", max_cells: int = 200_000):
        """Certify delta-density of accepted samples in the infinity norm."""
        from .coverage import check_delta_dense
        return check_delta_dense(self, points, delta, epsilon, norm=norm,
                                 method=method, max_cells=max_cells)

    def plot_projection(self, indices=None, *, basis=None, ranges=None,
                        plot_range=None, samples=15, max_sample_points=50_000,
                        max_candidates=None, seed=0, default_window=(-2, 2),
                        **kwds):
        """Return an approximate Sage 2D/3D projection graphic.

        If no finite window is supplied, simple coordinate bounds are inferred
        and remaining directions use a configurable default visualization
        window.  Sampling and plotting never certify that the set is empty.
        """
        from .plotting import plot_projection
        return plot_projection(self, indices=indices, basis=basis, ranges=ranges,
                               plot_range=plot_range, samples=samples,
                               max_sample_points=max_sample_points,
                               max_candidates=max_candidates, seed=seed,
                               default_window=default_window, **kwds)

