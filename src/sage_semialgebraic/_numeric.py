"""Small validation helpers shared by public numerical APIs."""

from __future__ import annotations

from numbers import Integral


def as_python_integer(value: object, name: str, *, minimum: int | None = None) -> int:
    """Validate a Python or Sage integer without truncating other numbers."""
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    result = int(value)
    if minimum is not None and result < minimum:
        if minimum == 1:
            raise ValueError(f"{name} must be positive")
        raise ValueError(f"{name} must be an integer at least {minimum}")
    return result

