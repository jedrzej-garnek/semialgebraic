"""Headless Sage graphics for explicitly approximate projections."""

from __future__ import annotations

from itertools import islice, product
from typing import Sequence

from sage.all import Graphics, QQ, RDF, matrix, point, point3d, text3d, vector
from sage.plot.plot3d.base import Graphics3dGroup


def _normalise_ranges(dimension: int, ranges, plot_range):
    if ranges is not None and plot_range is not None:
        raise ValueError("give only one of ranges and plot_range")
    value = ranges if ranges is not None else plot_range
    if value is None:
        raise ValueError("explicit finite ambient ranges are required")
    if (len(value) == 2 and not hasattr(value[0], "__len__")):
        value = [value] * dimension
    if len(value) != dimension:
        raise ValueError(f"expected {dimension} ambient plotting ranges")
    result = []
    for bounds in value:
        if len(bounds) != 2 or bounds[0] >= bounds[1]:
            raise ValueError("every plotting range must be a pair (lower, upper)")
        result.append((QQ(bounds[0]), QQ(bounds[1])))
    return tuple(result)


def _grid(ranges, samples: int, cap: int):
    axes = [tuple(lo + (hi - lo) * QQ(i) / (samples - 1)
                  for i in range(samples)) for lo, hi in ranges]
    return islice(product(*axes), cap)


def _sample_set(semiset, ranges, samples: int, cap: int):
    if semiset._kind == "basic":
        output = []
        for candidate in _grid(ranges, samples, cap):
            if semiset.contains(candidate):
                output.append(tuple(candidate))
        return output
    if semiset._kind == "product":
        left, right = semiset._operands
        n = left.ambient_dimension
        per_side = max(2, int(cap ** 0.5))
        left_points = _sample_set(left, ranges[:n], samples, per_side)
        right_points = _sample_set(right, ranges[n:], samples, per_side)
        return [a + b for a, b in islice(product(left_points, right_points), cap)]
    if semiset._kind == "join":
        left, right = semiset._operands
        per_side = max(2, int((cap / max(samples, 2)) ** 0.5))
        left_points = _sample_set(left, ranges, samples, per_side)
        right_points = _sample_set(right, ranges, samples, per_side)
        output = []
        for a, b in product(left_points, right_points):
            for i in range(samples):
                t = QQ(i) / (samples - 1)
                candidate = tuple(t * x + (1 - t) * y for x, y in zip(a, b))
                if all(lo <= x <= hi for x, (lo, hi) in zip(candidate, ranges)):
                    output.append(candidate)
                    if len(output) >= cap:
                        return output
        return output
    return []


def _projection(ambient_dimension: int, indices, basis):
    if (indices is None) == (basis is None):
        if indices is None:
            indices = tuple(range(min(ambient_dimension, 2)))
        else:
            raise ValueError("give exactly one of indices and basis")
    if basis is None:
        indices = tuple(int(i) for i in indices)
        if len(indices) not in (2, 3):
            raise ValueError("projection dimension must be 2 or 3")
        if len(set(indices)) != len(indices) or any(
                i < 0 or i >= ambient_dimension for i in indices):
            raise ValueError("projection indices must be distinct ambient coordinates")
        return len(indices), lambda p: tuple(RDF(p[i]) for i in indices)
    supplied = tuple(tuple(v) for v in basis)
    if len(supplied) not in (2, 3):
        raise ValueError("basis must contain two or three vectors")
    if any(len(v) != ambient_dimension for v in supplied):
        raise ValueError("each basis vector must have ambient dimension entries")
    raw = [vector(RDF, v) for v in supplied]
    if matrix(RDF, supplied).rank() != len(supplied):
        raise ValueError("basis vectors must be linearly independent")
    orthonormal = []
    for v in raw:
        w = v
        for q in orthonormal:
            w = w - w.dot_product(q) * q
        orthonormal.append(w / w.norm())
    return len(supplied), lambda p: tuple(
        vector(RDF, p).dot_product(q) for q in orthonormal)


def plot_projection(semiset, indices=None, *, basis=None, ranges=None,
                    plot_range=None, samples=15, max_sample_points=50_000,
                    point_size=8, color="steelblue", **kwds):
    """Draw a numerical sample of an orthogonal projection in 2D or 3D.

    The returned object is a Sage ``Graphics``/``Graphics3d`` object.  It is a
    visualization only, never evidence for an exact containment or coverage
    statement.  Finite ambient ranges are mandatory.
    """
    if not isinstance(samples, int) or samples < 2:
        raise ValueError("samples must be an integer at least 2")
    if not isinstance(max_sample_points, int) or max_sample_points <= 0:
        raise ValueError("max_sample_points must be positive")
    ambient_ranges = _normalise_ranges(semiset.ambient_dimension, ranges,
                                      plot_range)
    dimension, project = _projection(semiset.ambient_dimension, indices, basis)
    points = [project(p) for p in _sample_set(
        semiset, ambient_ranges, samples, max_sample_points)]
    if dimension == 2:
        if not points:
            graphic = Graphics()
        else:
            graphic = point(points, size=point_size, color=color,
                            legend_label="approximate projection", **kwds)
        graphic._extra_kwds["sage_semialgebraic_note"] = (
            "approximate sampled projection; not a proof")
        return graphic
    if not points:
        graphic3d = Graphics3dGroup([])
    else:
        graphic3d = point3d(points, size=point_size, color=color, **kwds)
        minima = tuple(min(p[i] for p in points) for i in range(3))
        graphic3d += text3d("approximation", minima, color="black")
    return graphic3d

