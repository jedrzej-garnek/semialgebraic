"""Visualization sampling with exact validation of returned points.

Sampling is never used as an emptiness proof.  Basic-set samples are checked
against exact constraints.  Polynomial hypersurfaces are sampled by exact
specialization and real-algebraic root isolation in Sage's ``AA`` after exact
elimination of linear equalities over ``QQ``.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import chain, product
from math import gcd, isqrt
from random import Random
from typing import Iterable, Iterator, Sequence

from sage.all import AA, QQ, PolynomialRing, matrix

from .coverage import Interval, _linear_constraint


@dataclass(frozen=True)
class SamplingResult:
    """Internal result carrying sample points and non-proof diagnostics."""

    points: tuple[tuple[object, ...], ...]
    diagnostics: tuple[str, ...] = ()


def _primes(count: int) -> tuple[int, ...]:
    primes: list[int] = []
    candidate = 2
    while len(primes) < count:
        if all(candidate % p for p in primes if p * p <= candidate):
            primes.append(candidate)
        candidate += 1
    return tuple(primes)


def _radical_inverse(index: int, base: int):
    numerator = 0
    denominator = 1
    while index:
        index, digit = divmod(index, base)
        numerator = numerator * base + digit
        denominator *= base
    return QQ(numerator) / denominator


def _unique_points(points: Iterable[Sequence[object]], limit: int):
    seen = set()
    output = []
    for point_ in points:
        point_tuple = tuple(point_)
        try:
            key = point_tuple
            if key in seen:
                continue
        except TypeError:
            key = tuple(repr(x) for x in point_tuple)
            if key in seen:
                continue
        seen.add(key)
        output.append(point_tuple)
        if len(output) >= limit:
            break
    return tuple(output)


def candidate_points(ranges: Sequence[tuple[object, object]], samples: int,
                     cap: int, seed: int) -> Iterator[tuple[object, ...]]:
    """Yield a full grid or a window-wide deterministic Halton sample."""
    if cap <= 0:
        return
    dimension = len(ranges)
    full_size = samples ** dimension
    if full_size <= cap:
        axes = [tuple(lo + (hi - lo) * QQ(i) / (samples - 1)
                      for i in range(samples)) for lo, hi in ranges]
        yield from product(*axes)
        return

    seen: set[tuple[object, ...]] = set()
    emitted = 0
    centre = tuple((lo + hi) / 2 for lo, hi in ranges)
    anchors = [centre]
    if dimension <= 12 and 2 ** dimension <= max(2, cap // 8):
        anchors.extend(tuple(bounds[bit] for bounds, bit in zip(ranges, bits))
                       for bits in product((0, 1), repeat=dimension))
    for anchor in anchors:
        if anchor not in seen:
            seen.add(anchor)
            yield anchor
            emitted += 1
            if emitted >= cap:
                return

    bases = _primes(dimension)
    index = 1 + (abs(seed) * 104729) % 1_000_003
    while emitted < cap:
        candidate = tuple(
            lo + (hi - lo) * _radical_inverse(index, base)
            for (lo, hi), base in zip(ranges, bases))
        index += 1
        if candidate in seen:
            continue
        seen.add(candidate)
        yield candidate
        emitted += 1


def _spread_indices(total: int, count: int, seed: int) -> Iterator[int]:
    """Select unique indices across a large implicit Cartesian product."""
    if total <= count:
        yield from range(total)
        return
    offset = (abs(seed) * 104729 + total // 2) % total
    stride = max(1, total // count) + 2 * abs(seed) + 1
    while gcd(stride, total) != 1:
        stride += 1
    current = offset
    for _ in range(count):
        yield current
        current = (current + stride) % total


def _within_ranges(point_: Sequence[object], ranges) -> bool:
    return all(lo <= value <= hi for value, (lo, hi) in zip(point_, ranges))


def _total_degree(polynomial) -> int:
    """Return total degree for Sage univariate and multivariate polynomials."""
    total_degree = getattr(polynomial, "total_degree", None)
    return int(total_degree() if total_degree is not None else polynomial.degree())


def _linear_parameterization(semiset, linear_equalities):
    """Return affine coordinate expressions and free columns, or ``None``."""
    n = semiset.ambient_dimension
    generators = tuple(semiset._ring.gens())
    zero = {g: 0 for g in generators}
    rows, rhs = [], []
    for polynomial in linear_equalities:
        rows.append([QQ(polynomial.derivative(g)) for g in generators])
        rhs.append(-QQ(polynomial.subs(zero)))
    if not rows:
        parameter_ring = PolynomialRing(QQ, n, names=tuple(f"u{i}" for i in range(n)))
        return tuple(parameter_ring.gens()), tuple(range(n)), parameter_ring
    augmented = matrix(QQ, rows).augment(matrix(QQ, len(rows), 1, rhs))
    echelon = augmented.echelon_form()
    for row in range(echelon.nrows()):
        if all(echelon[row, column] == 0 for column in range(n)) and \
                echelon[row, n] != 0:
            return None
    pivot_columns = []
    pivot_rows = {}
    for row in range(echelon.nrows()):
        for column in range(n):
            if echelon[row, column] != 0:
                pivot_columns.append(column)
                pivot_rows[column] = row
                break
    free_columns = tuple(i for i in range(n) if i not in pivot_columns)
    parameter_count = max(1, len(free_columns))
    parameter_ring = PolynomialRing(
        QQ, parameter_count, names=tuple(f"u{i}" for i in range(parameter_count)))
    parameters = tuple(parameter_ring.gens())[:len(free_columns)]
    expressions: list[object] = [parameter_ring.zero() for _ in range(n)]
    for index, column in enumerate(free_columns):
        expressions[column] = parameters[index]
    for pivot in pivot_columns:
        row = pivot_rows[pivot]
        coefficient = echelon[row, pivot]
        expression = parameter_ring(echelon[row, n] / coefficient)
        for index, column in enumerate(free_columns):
            expression -= (echelon[row, column] / coefficient) * parameters[index]
        expressions[pivot] = expression
    return tuple(expressions), free_columns, parameter_ring


def _evaluate_expressions(expressions, parameters, values):
    substitutions = dict(zip(parameters, values))
    return tuple(expression.subs(substitutions) for expression in expressions)


def _hypersurface_points(polynomial, ranges, samples: int, budget: int,
                         seed: int) -> Iterator[tuple[object, ...]]:
    """Yield exact points on one polynomial hypersurface within a box."""
    generators = tuple(polynomial.parent().gens())
    active = [i for i, generator in enumerate(generators)
              if polynomial.degree(generator) > 0]
    if not active:
        return
    per_variable = max(1, budget // len(active))
    for turn, solve_index in enumerate(active):
        other_indices = tuple(i for i in range(len(generators)) if i != solve_index)
        other_ranges = tuple(ranges[i] for i in other_indices)
        fixed_candidates = candidate_points(
            other_ranges, samples, per_variable, seed + 7919 * (turn + 1)) \
            if other_indices else iter(((),))
        checked = 0
        for fixed in fixed_candidates:
            if checked >= per_variable:
                break
            checked += 1
            substitutions = {generators[i]: value
                             for i, value in zip(other_indices, fixed)}
            specialized = polynomial.subs(substitutions)
            if specialized == 0:
                lo, hi = ranges[solve_index]
                roots = (lo, (lo + hi) / 2, hi)
            elif not specialized.variables():
                continue
            else:
                try:
                    roots = specialized.univariate_polynomial().roots(
                        ring=AA, multiplicities=False)
                except (ArithmeticError, TypeError, ValueError, NotImplementedError):
                    continue
            for root in roots:
                if not ranges[solve_index][0] <= root <= ranges[solve_index][1]:
                    continue
                point_ = [None] * len(generators)
                point_[solve_index] = root
                for i, value in zip(other_indices, fixed):
                    point_[i] = value
                yield tuple(point_)


def _sample_basic(semiset, ranges, samples: int, candidate_cap: int,
                  point_cap: int, seed: int) -> SamplingResult:
    equalities = tuple(p for p in semiset._constraints["eq"] if p != 0)
    if not equalities:
        accepted = (candidate for candidate in candidate_points(
            ranges, samples, candidate_cap, seed) if semiset.contains(candidate))
        return SamplingResult(_unique_points(accepted, point_cap))

    linear = tuple(p for p in equalities if _total_degree(p) <= 1)
    nonlinear = tuple(p for p in equalities if _total_degree(p) > 1)
    parameterization = _linear_parameterization(semiset, linear)
    if parameterization is None:
        return SamplingResult((), ("the linear equality subsystem is inconsistent",))
    expressions, free_columns, parameter_ring = parameterization
    parameters = tuple(parameter_ring.gens())[:len(free_columns)]
    parameter_ranges = tuple(ranges[i] for i in free_columns)

    def lift(values):
        point_ = _evaluate_expressions(expressions, parameters, values)
        if _within_ranges(point_, ranges) and semiset.contains(point_):
            return point_
        return None

    if not free_columns:
        point_ = lift(())
        return SamplingResult(() if point_ is None else (point_,))

    transformed = []
    for polynomial in nonlinear:
        image = parameter_ring(polynomial(*expressions))
        if image == 0:
            continue
        if not image.variables():
            return SamplingResult((), ("the equality subsystem is inconsistent",))
        transformed.append(image)

    if not transformed:
        centre = tuple((lo + hi) / 2 for lo, hi in parameter_ranges)
        parameter_candidates = chain(
            (centre,), candidate_points(
                parameter_ranges, samples, candidate_cap - 1, seed))
        lifted = (lift(values) for values in parameter_candidates)
        return SamplingResult(_unique_points(
            (point_ for point_ in lifted if point_ is not None), point_cap))

    per_equation = max(1, candidate_cap // len(transformed))

    def exact_points():
        for index, driver in enumerate(transformed):
            for values in _hypersurface_points(
                    driver, parameter_ranges, samples, per_equation,
                    seed + 65537 * (index + 1)):
                point_ = lift(values)
                if point_ is not None:
                    yield point_

    points = _unique_points(exact_points(), point_cap)
    diagnostics = ()
    if not points:
        diagnostics = (
            "exact specialization found no sample points; this is not an emptiness proof",
        )
    elif len(transformed) > 1:
        diagnostics = (
            "multiple nonlinear equalities were sampled by hypersurface specialization; "
            "coverage of every component is not guaranteed",
        )
    return SamplingResult(points, diagnostics)


def _join_leaves(semiset):
    """Flatten nested joins, retaining their exact set-theoretic meaning.

    A join of k sets is the set of convex combinations containing one point
    from each of the k sets. This does NOT convexify any individual operand.
    """
    if semiset._kind == "join":
        left, right = semiset._operands
        yield from _join_leaves(left)
        yield from _join_leaves(right)
    else:
        yield semiset


def _simplex_weights(count: int, rng: Random):
    """Return exact rational weights, approximately uniform on the simplex."""
    if count == 1:
        return (QQ.one(),)
    denominator = 1 << 32
    cuts = (0, *sorted(rng.sample(range(1, denominator), count - 1)),
            denominator)
    return tuple(QQ(cuts[i + 1] - cuts[i]) / denominator
                 for i in range(count))


def _sample_join(semiset, ranges, samples: int, candidate_cap: int,
                 point_cap: int, seed: int) -> SamplingResult:
    """Sample a flattened join without creating artificial intermediate fans.

    Samples from different operands and simplex weights are selected
    independently, which avoids the correlated-index problem of
    ``_spread_indices`` on a flattened ``(left, right, t)`` product.
    """
    leaves = tuple(_join_leaves(semiset))
    leaf_count = len(leaves)
    leaf_point_cap = min(point_cap, max(samples, 4 * isqrt(point_cap)))
    leaf_candidate_cap = max(1, candidate_cap // leaf_count)

    results = tuple(
        sample_set(leaf, ranges, samples, leaf_candidate_cap, leaf_point_cap,
                   seed + 7919 * (i + 1))
        for i, leaf in enumerate(leaves)
    )
    diagnostics = tuple(note for result in results
                        for note in result.diagnostics)
    if any(not result.points for result in results):
        return SamplingResult(
            (), diagnostics + (
                "a join operand had no sampled points; "
                "this is not an emptiness proof",
            ))

    operands = tuple(result.points for result in results)
    rng = Random(seed)
    denominator = 1 << 32
    dimension = semiset.ambient_dimension

    def combinations():
        # Include a small sample of the individual operands (simplex vertices).
        # Reserve most of the display budget for mixtures, even when a leaf has
        # very many sample points.
        from itertools import islice

        def operand_anchors():
            for index in range(max(map(len, operands))):
                for operand_points in operands:
                    if index < len(operand_points):
                        yield operand_points[index]

        yield from islice(operand_anchors(), max(0, point_cap // 10))

        for attempt in range(max(candidate_cap, point_cap)):
            chosen = tuple(rng.choice(operand_points)
                           for operand_points in operands)
            if attempt % 4 == 0:
                # Also sample pairwise joins: helpful for displaying the
                # boundary of cones and related joins with point operands.
                first, second = rng.sample(range(leaf_count), 2)
                t = QQ(rng.randrange(1, denominator)) / denominator
                weights = [QQ.zero()] * leaf_count
                weights[first] = t
                weights[second] = 1 - t
            else:
                weights = _simplex_weights(leaf_count, rng)
            candidate = tuple(
                sum(weights[i] * chosen[i][coord]
                    for i in range(leaf_count))
                for coord in range(dimension)
            )
            if _within_ranges(candidate, ranges):
                yield candidate

    return SamplingResult(_unique_points(combinations(), point_cap), diagnostics)


def sample_set(semiset, ranges, samples: int, candidate_cap: int,
               point_cap: int, seed: int = 0) -> SamplingResult:
    """Sample a set recursively while preserving exact join/product semantics."""
    if semiset._kind == "basic":
        return _sample_basic(semiset, ranges, samples, candidate_cap, point_cap, seed)

    if semiset._kind == "product":
        left, right = semiset._operands
        n = left.ambient_dimension
        side_limit = min(point_cap, max(samples, 2 * isqrt(point_cap) + 1))
        left_result = sample_set(left, ranges[:n], samples,
                                 max(1, candidate_cap // 2), side_limit, seed + 1)
        right_result = sample_set(right, ranges[n:], samples,
                                  max(1, candidate_cap // 2), side_limit, seed + 2)
        left_points, right_points = left_result.points, right_result.points
        total = len(left_points) * len(right_points)

        def combinations():
            if not right_points:
                return
            for flat in _spread_indices(total, point_cap, seed + 3):
                i, j = divmod(flat, len(right_points))
                yield left_points[i] + right_points[j]

        return SamplingResult(
            _unique_points(combinations(), point_cap),
            left_result.diagnostics + right_result.diagnostics)

    if semiset._kind == "join":
        return _sample_join(semiset, ranges, samples, candidate_cap,
                            point_cap, seed)

    return SamplingResult((), ("unsupported symbolic set node for plotting",))


def _basic_coordinate_bounds(semiset) -> tuple[Interval, ...]:
    intervals = [Interval() for _ in range(semiset.ambient_dimension)]
    for relation in ("eq", "ge", "gt", "le", "lt"):
        for polynomial in semiset._constraints[relation]:
            parsed = _linear_constraint(polynomial)
            if parsed is False:
                continue
            index, coefficient, constant = parsed
            if index is None or coefficient == 0:
                continue
            boundary = -constant / coefficient
            effective = relation
            if coefficient < 0:
                effective = {"ge": "le", "gt": "lt", "le": "ge",
                             "lt": "gt", "eq": "eq"}[relation]
            current = intervals[index]
            lower, upper = current.lower, current.upper
            lower_closed, upper_closed = current.lower_closed, current.upper_closed
            if effective == "eq":
                lower = upper = boundary
                lower_closed = upper_closed = True
            elif effective in ("ge", "gt") and (lower is None or boundary > lower):
                lower, lower_closed = boundary, effective == "ge"
            elif effective in ("le", "lt") and (upper is None or boundary < upper):
                upper, upper_closed = boundary, effective == "le"
            intervals[index] = Interval(lower, upper, lower_closed, upper_closed)
    return tuple(intervals)


def infer_ranges(semiset, default_window=(-2, 2)):
    """Infer a finite visualization window and explanatory diagnostics."""
    default_lower, default_upper = (QQ(default_window[0]), QQ(default_window[1]))
    if default_lower >= default_upper:
        raise ValueError("default_window must be a pair (lower, upper)")
    width = default_upper - default_lower
    if semiset._kind == "basic":
        intervals = _basic_coordinate_bounds(semiset)
    elif semiset._kind == "product":
        left_ranges, left_notes = infer_ranges(semiset._operands[0], default_window)
        right_ranges, right_notes = infer_ranges(semiset._operands[1], default_window)
        return left_ranges + right_ranges, left_notes + right_notes
    elif semiset._kind == "join":
        left_ranges, left_notes = infer_ranges(semiset._operands[0], default_window)
        right_ranges, right_notes = infer_ranges(semiset._operands[1], default_window)
        combined = tuple((min(a[0], b[0]), max(a[1], b[1]))
                         for a, b in zip(left_ranges, right_ranges))
        return combined, left_notes + right_notes
    else:
        intervals = tuple(Interval() for _ in range(semiset.ambient_dimension))

    ranges = []
    notes = []
    for coordinate, interval in enumerate(intervals):
        if interval.lower is not None and interval.upper is not None:
            if interval.lower < interval.upper:
                ranges.append((interval.lower, interval.upper))
            else:
                padding = width / 4
                ranges.append((interval.lower - padding, interval.upper + padding))
        elif interval.lower is not None:
            ranges.append((interval.lower, interval.lower + width))
            notes.append(f"coordinate {coordinate} is unbounded above; a finite window was used")
        elif interval.upper is not None:
            ranges.append((interval.upper - width, interval.upper))
            notes.append(f"coordinate {coordinate} is unbounded below; a finite window was used")
        else:
            ranges.append((default_lower, default_upper))
            notes.append(f"coordinate {coordinate} uses the default visualization window")
    return tuple(ranges), tuple(notes)
