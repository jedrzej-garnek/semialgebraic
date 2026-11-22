"""Certified elementary coverage algorithms.

The box backend partitions every coordinate at all sample-box and domain
endpoints.  Open intervals and singleton boundary cells are both represented,
which is essential: testing only a lattice of cube origins is not a proof.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import shutil
from typing import Iterable, Sequence

from sage.all import QQ

from ._numeric import as_python_integer
from .result import CoverageResult, CoverageStatus


@dataclass(frozen=True)
class Interval:
    lower: object | None = None
    upper: object | None = None
    lower_closed: bool = False
    upper_closed: bool = False

    def contains(self, value: object) -> bool:
        if self.lower is not None:
            if value < self.lower or (value == self.lower and not self.lower_closed):
                return False
        if self.upper is not None:
            if value > self.upper or (value == self.upper and not self.upper_closed):
                return False
        return True

    @property
    def empty(self) -> bool:
        if self.lower is None or self.upper is None:
            return False
        return (self.lower > self.upper or
                (self.lower == self.upper and
                 not (self.lower_closed and self.upper_closed)))


@dataclass(frozen=True)
class AxisBox:
    intervals: tuple[Interval, ...]
    empty: bool = False


def _constant_relation(value: object, relation: str) -> bool:
    return {
        "eq": value == 0,
        "ge": value >= 0,
        "gt": value > 0,
        "le": value <= 0,
        "lt": value < 0,
    }[relation]


def _tighten_lower(interval: Interval, value: object, closed: bool) -> Interval:
    if interval.lower is None or value > interval.lower:
        return Interval(value, interval.upper, closed, interval.upper_closed)
    if value == interval.lower:
        return Interval(value, interval.upper,
                        interval.lower_closed and closed, interval.upper_closed)
    return interval


def _tighten_upper(interval: Interval, value: object, closed: bool) -> Interval:
    if interval.upper is None or value < interval.upper:
        return Interval(interval.lower, value, interval.lower_closed, closed)
    if value == interval.upper:
        return Interval(interval.lower, value, interval.lower_closed,
                        interval.upper_closed and closed)
    return interval


def _linear_constraint(polynomial: object):
    """Return ``(variable index or None, a, b)`` for ``a*x+b``."""
    ring = polynomial.parent()
    variables = tuple(polynomial.variables())
    zero_subs = {g: 0 for g in ring.gens()}
    if not variables:
        return None, QQ(0), QQ(polynomial)
    if len(variables) != 1:
        return False
    variable = variables[0]
    derivative = polynomial.derivative(variable)
    if tuple(derivative.variables()):
        return False
    a = QQ(derivative)
    b = QQ(polynomial.subs(zero_subs))
    if polynomial != a * variable + b:
        return False
    index = tuple(ring.gens()).index(variable)
    return index, a, b


def extract_axis_box(semiset) -> AxisBox | None:
    """Recognize an exact Cartesian product of rational intervals."""
    if semiset._kind == "product":
        left = extract_axis_box(semiset._operands[0])
        right = extract_axis_box(semiset._operands[1])
        if left is None or right is None:
            return None
        return AxisBox(left.intervals + right.intervals, left.empty or right.empty)
    if semiset._kind == "join":
        a = singleton_point(semiset._operands[0])
        b = singleton_point(semiset._operands[1])
        if a is not None and b is not None:
            if semiset.ambient_dimension == 1:
                lo, hi = sorted((QQ(a[0]), QQ(b[0])))
                return AxisBox((Interval(lo, hi, True, True),))
            if a == b:
                return AxisBox(tuple(Interval(QQ(x), QQ(x), True, True)
                                         for x in a))
        return None
    if semiset._kind != "basic":
        return None
    intervals = [Interval() for _ in range(semiset.ambient_dimension)]
    for relation in ("eq", "ge", "gt", "le", "lt"):
        for polynomial in semiset._constraints[relation]:
            parsed = _linear_constraint(polynomial)
            if parsed is False:
                return None
            index, a, b = parsed
            if index is None or a == 0:
                if not _constant_relation(b, relation):
                    return AxisBox(tuple(intervals), empty=True)
                continue
            boundary = -b / a
            effective = relation
            if a < 0:
                effective = {"ge": "le", "gt": "lt", "le": "ge",
                             "lt": "gt", "eq": "eq"}[relation]
            current = intervals[index]
            if effective == "eq":
                current = _tighten_lower(current, boundary, True)
                current = _tighten_upper(current, boundary, True)
            elif effective in ("ge", "gt"):
                current = _tighten_lower(current, boundary, effective == "ge")
            else:
                current = _tighten_upper(current, boundary, effective == "le")
            intervals[index] = current
            if current.empty:
                return AxisBox(tuple(intervals), empty=True)
    return AxisBox(tuple(intervals), any(i.empty for i in intervals))


def singleton_point(semiset) -> tuple[object, ...] | None:
    box = extract_axis_box(semiset)
    if box is None or box.empty:
        return None
    result = []
    for interval in box.intervals:
        if (interval.lower is None or interval.upper is None or
                interval.lower != interval.upper or
                not interval.lower_closed or not interval.upper_closed):
            return None
        result.append(interval.lower)
    return tuple(result)


def _has_nonzero_equality(semiset) -> bool:
    if semiset._kind == "basic":
        return any(polynomial != 0 for polynomial in semiset._constraints["eq"])
    if semiset._kind == "product":
        return any(_has_nonzero_equality(s) for s in semiset._operands)
    if semiset._kind == "join":
        a, b = (singleton_point(s) for s in semiset._operands)
        return (a is not None and b is not None and
                (a == b or semiset.ambient_dimension > 1))
    return False


def qepcad_available() -> bool:
    """Whether Sage's interface and the external QEPCAD executable both exist."""
    if shutil.which("qepcad") is None:
        return False
    try:
        import sage.interfaces.qepcad  # noqa: F401
    except ImportError:
        return False
    return True


def _positive_rational(value: object, name: str):
    try:
        result = QQ(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise TypeError(f"{name} must be an exact real value coercible to QQ") from exc
    if result <= 0:
        raise ValueError(f"{name} must be positive")
    return result


def _validate_epsilon(epsilon: object) -> None:
    try:
        if epsilon < 0:
            raise ValueError("epsilon must be nonnegative")
    except TypeError as exc:
        raise TypeError("epsilon must be an ordered real number") from exc


def _materialize_points(semiset, points: Iterable[Sequence[object]]):
    materialized = []
    for point in points:
        materialized.append(semiset._validate_point(point))
    return tuple(materialized)


def _accepted_rational_points(semiset, points, epsilon):
    accepted, rejected = [], []
    for point in points:
        try:
            inside = semiset.contains(point, epsilon=epsilon)
        except NotImplementedError:
            return None, None
        if not inside:
            rejected.append(point)
            continue
        try:
            accepted.append(tuple(QQ(x) for x in point))
        except (TypeError, ValueError, OverflowError):
            return None, None
    return tuple(accepted), tuple(rejected)


def _cube_origin_box(box: AxisBox, delta: object) -> AxisBox:
    intervals = []
    for interval in box.intervals:
        upper = None if interval.upper is None else interval.upper - delta
        intervals.append(Interval(interval.lower, upper,
                                  interval.lower_closed,
                                  interval.upper_closed))
    return AxisBox(tuple(intervals), box.empty or any(i.empty for i in intervals))


def _partition_representatives(domain: Interval,
                               cover_intervals: Sequence[Interval]):
    endpoints = set()
    if domain.lower is not None:
        endpoints.add(domain.lower)
    if domain.upper is not None:
        endpoints.add(domain.upper)
    for interval in cover_intervals:
        if interval.lower is not None:
            endpoints.add(interval.lower)
        if interval.upper is not None:
            endpoints.add(interval.upper)
    ordered = sorted(endpoints)
    candidates = []
    if not ordered:
        candidates.append(QQ(0))
    else:
        candidates.append(ordered[0] - 1)
        for i, endpoint in enumerate(ordered):
            candidates.append(endpoint)
            if i + 1 < len(ordered):
                candidates.append((endpoint + ordered[i + 1]) / 2)
        candidates.append(ordered[-1] + 1)
    return tuple(value for value in candidates if domain.contains(value))


def _cover_box(domain: AxisBox, cover_boxes: Sequence[AxisBox], *,
               max_cells: int, backend: str, claim: str) -> CoverageResult:
    if max_cells <= 0:
        raise ValueError("max_cells must be positive")
    if domain.empty:
        return CoverageResult(CoverageStatus.PROVED, backend,
                              f"{claim} is vacuous because the exact domain is empty",
                              metadata={"cells_checked": 0})
    per_coordinate = []
    for i, interval in enumerate(domain.intervals):
        reps = _partition_representatives(
            interval, [box.intervals[i] for box in cover_boxes])
        if not reps:
            return CoverageResult(CoverageStatus.PROVED, backend,
                                  f"{claim} is vacuous because the domain is empty",
                                  metadata={"cells_checked": 0})
        per_coordinate.append(reps)
    cell_count = 1
    for reps in per_coordinate:
        cell_count *= len(reps)
    if cell_count > max_cells:
        return CoverageResult(
            CoverageStatus.UNKNOWN, backend,
            f"exact endpoint decomposition needs {cell_count} cells, exceeding "
            f"max_cells={max_cells}",
            metadata={"required_cells": cell_count, "max_cells": max_cells})
    checked = 0
    for candidate in product(*per_coordinate):
        checked += 1
        covered = any(all(interval.contains(value)
                          for interval, value in zip(box.intervals, candidate))
                      for box in cover_boxes)
        if not covered:
            return CoverageResult(
                CoverageStatus.DISPROVED, backend,
                f"an exact uncovered {claim} witness was found",
                witness=tuple(candidate),
                metadata={"cells_checked": checked, "witness_verified": True})
    return CoverageResult(
        CoverageStatus.PROVED, backend,
        f"all {cell_count} endpoint-decomposition cells are covered exactly",
        metadata={"cells_checked": cell_count})


def _unsupported(method: str, subject: str) -> CoverageResult:
    if method == "qepcad" and not qepcad_available():
        reason = "Sage's QEPCAD module is importable but the qepcad executable is absent"
    elif qepcad_available():
        reason = (f"{subject} is outside the certified elementary backend; the "
                  "QEPCAD executable was detected but formula translation is not implemented")
    else:
        reason = (f"{subject} is outside the certified rational-box backend and no "
                  "QEPCAD executable is available")
    return CoverageResult(CoverageStatus.UNKNOWN, "none", reason)


def check_cube_hitting(semiset, points, delta, epsilon=0, *, method="auto",
                       max_cells=200_000) -> CoverageResult:
    delta = _positive_rational(delta, "delta")
    max_cells = as_python_integer(max_cells, "max_cells", minimum=1)
    _validate_epsilon(epsilon)
    if method not in ("auto", "box", "qepcad"):
        raise ValueError("method must be 'auto', 'box', or 'qepcad'")
    materialized = _materialize_points(semiset, points)
    if _has_nonzero_equality(semiset):
        return CoverageResult(
            CoverageStatus.PROVED, "empty-interior-certificate",
            "a nonzero polynomial equality forces empty ambient interior, so no "
            "positive-side-length ambient cube is contained in the set",
            metadata={"vacuous": True})
    box = extract_axis_box(semiset)
    if box is None or method == "qepcad":
        return _unsupported(method, "cube hitting for this formula")
    accepted, rejected = _accepted_rational_points(semiset, materialized, epsilon)
    if accepted is None:
        return _unsupported(method, "validation of these sample coordinates")
    origins = _cube_origin_box(box, delta)
    hit_boxes = tuple(AxisBox(tuple(
        Interval(p - delta, p, True, True) for p in point))
        for point in accepted)
    result = _cover_box(origins, hit_boxes, max_cells=max_cells,
                        backend="rational-axis-box", claim="cube-origin")
    metadata = dict(result.metadata)
    metadata.update({"accepted_points": len(accepted),
                     "rejected_points": len(rejected)})
    return CoverageResult(result.status, result.backend, result.reason,
                          result.witness, metadata)


def check_delta_dense(semiset, points, delta, epsilon=0, *, norm="linf",
                      method="auto", max_cells=200_000) -> CoverageResult:
    delta = _positive_rational(delta, "delta")
    max_cells = as_python_integer(max_cells, "max_cells", minimum=1)
    _validate_epsilon(epsilon)
    if norm != "linf":
        raise NotImplementedError("only the certified infinity-norm backend is implemented")
    if method not in ("auto", "box", "qepcad"):
        raise ValueError("method must be 'auto', 'box', or 'qepcad'")
    materialized = _materialize_points(semiset, points)
    box = extract_axis_box(semiset)
    if box is None or method == "qepcad":
        return _unsupported(method, "delta-density for this formula")
    accepted, rejected = _accepted_rational_points(semiset, materialized, epsilon)
    if accepted is None:
        return _unsupported(method, "validation of these sample coordinates")
    balls = tuple(AxisBox(tuple(
        Interval(p - delta, p + delta, True, True) for p in point))
        for point in accepted)
    result = _cover_box(box, balls, max_cells=max_cells,
                        backend="rational-axis-box", claim="set-point")
    metadata = dict(result.metadata)
    metadata.update({"accepted_points": len(accepted),
                     "rejected_points": len(rejected)})
    return CoverageResult(result.status, result.backend, result.reason,
                          result.witness, metadata)

