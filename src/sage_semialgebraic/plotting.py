"""Headless Sage graphics for explicitly approximate projections."""

from __future__ import annotations

import warnings

from sage.all import Graphics, QQ, RDF, matrix, point, point3d, text3d, vector
from sage.plot.plot3d.base import Graphics3dGroup

from ._numeric import as_python_integer
from .sampling import infer_ranges, sample_set


class PlottingWindowWarning(UserWarning):
    """A finite visualization window was inferred or substituted."""


class SamplingWarning(UserWarning):
    """The visualization sampler reports an incomplete strategy or no points."""


def _normalise_ranges(dimension: int, ranges, plot_range):
    if ranges is not None and plot_range is not None:
        raise ValueError("give only one of ranges and plot_range")
    value = ranges if ranges is not None else plot_range
    if value is None:
        return None
    if len(value) == 2 and not hasattr(value[0], "__len__"):
        value = [value] * dimension
    if len(value) != dimension:
        raise ValueError(f"expected {dimension} ambient plotting ranges")
    result = []
    for bounds in value:
        if len(bounds) != 2 or bounds[0] >= bounds[1]:
            raise ValueError("every plotting range must be a pair (lower, upper)")
        result.append((QQ(bounds[0]), QQ(bounds[1])))
    return tuple(result)


def _projection(ambient_dimension: int, indices, basis):
    if (indices is None) == (basis is None):
        if indices is None:
            indices = tuple(range(min(ambient_dimension, 2)))
        else:
            raise ValueError("give exactly one of indices and basis")
    if basis is None:
        indices = tuple(as_python_integer(i, "projection index", minimum=0)
                        for i in indices)
        if len(indices) not in (2, 3):
            raise ValueError("projection dimension must be 2 or 3")
        if len(set(indices)) != len(indices) or any(
                i >= ambient_dimension for i in indices):
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
                    max_candidates=None, seed=0, default_window=(-2, 2),
                    point_size=8, color="steelblue", **kwds):
    """Draw a reproducible numerical sample of an orthogonal projection.

    This returns a Sage 2D or 3D graphics object and never calls ``show()``.
    The drawing is an approximation, not an emptiness, membership, or coverage
    proof.  Without explicit ranges, coordinate bounds are inferred when
    possible and otherwise ``default_window`` is used; warnings describe
    clipped or defaulted directions.

    ``max_candidates`` bounds sampling work and ``max_sample_points`` bounds
    displayed points.  Python integers and Sage ``Integer`` values are valid.
    """
    samples = as_python_integer(samples, "samples", minimum=2)
    max_sample_points = as_python_integer(
        max_sample_points, "max_sample_points", minimum=1)
    if max_candidates is None:
        max_candidates = max_sample_points
    max_candidates = as_python_integer(max_candidates, "max_candidates", minimum=1)
    seed = as_python_integer(seed, "seed")
    ambient_ranges = _normalise_ranges(semiset.ambient_dimension, ranges, plot_range)
    if ambient_ranges is None:
        ambient_ranges, notes = infer_ranges(semiset, default_window)
        if notes:
            warnings.warn(
                "automatic plotting window: " + "; ".join(dict.fromkeys(notes)),
                PlottingWindowWarning, stacklevel=2)
    dimension, project = _projection(semiset.ambient_dimension, indices, basis)
    sampled = sample_set(semiset, ambient_ranges, samples, max_candidates,
                         max_sample_points, seed)
    for diagnostic in dict.fromkeys(sampled.diagnostics):
        warnings.warn(diagnostic, SamplingWarning, stacklevel=2)
    points = [project(p) for p in sampled.points]
    if not points:
        warnings.warn(
            "the visualization sampler found no points in the plotting window; "
            "this is not an emptiness proof", SamplingWarning, stacklevel=2)
    if dimension == 2:
        if not points:
            return Graphics()
        # legend_label is renderer-supported.  Unknown keys in _extra_kwds are
        # forwarded by Sage to matplotlib and must not be used as metadata.
        return point(points, size=point_size, color=color,
                     legend_label="approximate projection", **kwds)
    if not points:
        return Graphics3dGroup([])
    graphic3d = point3d(points, size=point_size, color=color, **kwds)
    minima = tuple(min(p[i] for p in points) for i in range(3))
    graphic3d += text3d("approximation", minima, color="black")
    return graphic3d

