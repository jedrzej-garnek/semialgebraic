"""Immutable first-order formula nodes used by :mod:`sage_semialgebraic`.

Atoms always compare a Sage polynomial with zero.  Quantified variables are
stored by name; all construction in the public API uses fresh names, so a
formula can be moved to a larger polynomial ring without variable capture.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Iterator, Mapping


class Relation(str, Enum):
    EQ = "="
    GE = ">="
    GT = ">"
    LE = "<="
    LT = "<"


class Formula:
    """Marker base class for immutable formula nodes."""


@dataclass(frozen=True)
class Truth(Formula):
    value: bool


@dataclass(frozen=True)
class Atom(Formula):
    polynomial: object
    relation: Relation


@dataclass(frozen=True)
class And(Formula):
    children: tuple[Formula, ...]

    def __init__(self, children: Iterable[Formula]):
        flattened: list[Formula] = []
        for child in children:
            if isinstance(child, Truth):
                if not child.value:
                    flattened = [child]
                    break
                continue
            if isinstance(child, And):
                flattened.extend(child.children)
            else:
                flattened.append(child)
        object.__setattr__(self, "children", tuple(flattened))


@dataclass(frozen=True)
class Or(Formula):
    children: tuple[Formula, ...]

    def __init__(self, children: Iterable[Formula]):
        object.__setattr__(self, "children", tuple(children))


@dataclass(frozen=True)
class Not(Formula):
    child: Formula


@dataclass(frozen=True)
class Exists(Formula):
    variables: tuple[str, ...]
    body: Formula

    def __init__(self, variables: Iterable[str], body: Formula):
        object.__setattr__(self, "variables", tuple(str(v) for v in variables))
        object.__setattr__(self, "body", body)


def iter_atoms(formula: Formula) -> Iterator[Atom]:
    """Yield all atoms in depth-first order."""
    if isinstance(formula, Atom):
        yield formula
    elif isinstance(formula, (And, Or)):
        for child in formula.children:
            yield from iter_atoms(child)
    elif isinstance(formula, Not):
        yield from iter_atoms(formula.child)
    elif isinstance(formula, Exists):
        yield from iter_atoms(formula.body)


def polynomial_variable_names(formula: Formula) -> tuple[str, ...]:
    """Return generator names appearing in atom polynomial parents."""
    names: list[str] = []
    for atom in iter_atoms(formula):
        for generator in atom.polynomial.parent().gens():
            name = str(generator)
            if name not in names:
                names.append(name)
    return tuple(names)


def remap_polynomial(polynomial: object, target_ring: object,
                     mapping: Mapping[str, object]) -> object:
    """Copy ``polynomial`` exactly to ``target_ring`` using generator names."""
    source_gens = tuple(polynomial.parent().gens())
    images = tuple(mapping[str(g)] for g in source_gens)
    result = target_ring.zero()
    for exponent, coefficient in polynomial.dict().items():
        # Sage's multivariate ETuple deliberately is not a built-in tuple.
        # Univariate dictionaries, in contrast, use an integer exponent.
        try:
            exponent = tuple(exponent)
        except TypeError:
            exponent = (exponent,)
        monomial = target_ring(coefficient)
        for image, power in zip(images, exponent):
            monomial *= image ** power
        result += monomial
    return result


def remap_formula(formula: Formula, target_ring: object,
                  mapping: Mapping[str, object]) -> Formula:
    """Copy a formula into ``target_ring`` without changing its structure."""
    if isinstance(formula, Truth):
        return formula
    if isinstance(formula, Atom):
        return Atom(remap_polynomial(formula.polynomial, target_ring, mapping),
                    formula.relation)
    if isinstance(formula, And):
        return And(remap_formula(c, target_ring, mapping) for c in formula.children)
    if isinstance(formula, Or):
        return Or(remap_formula(c, target_ring, mapping) for c in formula.children)
    if isinstance(formula, Not):
        return Not(remap_formula(formula.child, target_ring, mapping))
    if isinstance(formula, Exists):
        return Exists((str(mapping[v]) for v in formula.variables),
                      remap_formula(formula.body, target_ring, mapping))
    raise TypeError(f"unknown formula node {type(formula)!r}")


def format_formula(formula: Formula) -> str:
    """Return a compact, deterministic mathematical representation."""
    if isinstance(formula, Truth):
        return "true" if formula.value else "false"
    if isinstance(formula, Atom):
        return f"({formula.polynomial} {formula.relation.value} 0)"
    if isinstance(formula, And):
        return "true" if not formula.children else "(" + " and ".join(
            format_formula(c) for c in formula.children) + ")"
    if isinstance(formula, Or):
        return "false" if not formula.children else "(" + " or ".join(
            format_formula(c) for c in formula.children) + ")"
    if isinstance(formula, Not):
        return f"not {format_formula(formula.child)}"
    if isinstance(formula, Exists):
        return f"exists {', '.join(formula.variables)}: {format_formula(formula.body)}"
    raise TypeError(type(formula))

