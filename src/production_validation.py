"""Numerical and structural validation helpers for production results."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .dispersion import D_star_Lambda


@dataclass(frozen=True)
class ConjugateSummary:
    real_count: int
    pair_count: int
    unmatched: tuple[complex, ...]


def conjugate_summary(
    roots: Iterable[complex], *, relative_tolerance: float = 1e-8, absolute_tolerance: float = 1e-10
) -> ConjugateSummary:
    pending = [complex(root) for root in roots]
    real_count = 0
    pairs = 0
    unmatched: list[complex] = []
    while pending:
        root = pending.pop(0)
        scale = max(1.0, abs(root))
        if abs(root.imag) <= absolute_tolerance + relative_tolerance * scale:
            real_count += 1
            continue
        target = root.conjugate()
        match = next((i for i, other in enumerate(pending)
            if abs(other - target) <= absolute_tolerance + relative_tolerance * max(1.0, abs(root), abs(other))), None)
        if match is None:
            unmatched.append(root)
        else:
            pending.pop(match)
            pairs += 1
    return ConjugateSummary(real_count, pairs, tuple(unmatched))


def verify_residuals(
    roots: Iterable[complex], Arho: float, Amu: float, AG: float, Lambda: float,
    *, relative_tolerance: float = 1e-8, absolute_tolerance: float = 1e-10,
) -> tuple[float, ...]:
    residuals = tuple(abs(complex(D_star_Lambda(root, Arho, Amu, AG, Lambda))) for root in roots)
    threshold = absolute_tolerance + relative_tolerance
    if any(not np.isfinite(value) or value > threshold for value in residuals):
        raise ValueError("one or more roots fail original-equation residual verification")
    return residuals


def spectra_match(
    first: Iterable[complex], second: Iterable[complex], *, relative_tolerance: float = 1e-7, absolute_tolerance: float = 1e-9
) -> bool:
    pending = [complex(value) for value in second]
    for root in first:
        root = complex(root)
        match = next((index for index, other in enumerate(pending)
            if abs(root - other) <= absolute_tolerance + relative_tolerance * max(1.0, abs(root), abs(other))), None)
        if match is None:
            return False
        pending.pop(match)
    return not pending
