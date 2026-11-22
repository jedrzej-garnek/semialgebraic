"""Result types for certified coverage questions."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType
from typing import Any, Mapping


class CoverageStatus(str, Enum):
    """Three-valued status of an exact coverage claim."""

    PROVED = "proved"
    DISPROVED = "disproved"
    UNKNOWN = "unknown"


@dataclass(frozen=True)
class CoverageResult:
    """Explainable result of a cube-hitting or delta-density query."""

    status: CoverageStatus
    backend: str
    reason: str
    witness: tuple[object, ...] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", MappingProxyType(dict(self.metadata)))

    @property
    def proved(self) -> bool:
        return self.status is CoverageStatus.PROVED

    @property
    def used_backend(self) -> str:
        """Descriptive backend name (alias retained for readable diagnostics)."""
        return self.backend

    @property
    def counterexample(self) -> tuple[object, ...] | None:
        """Exact disproving witness, when one is available."""
        return self.witness

