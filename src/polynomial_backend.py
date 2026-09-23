"""Selectable numerical backends for the unchanged canonical P14 polynomial."""

from __future__ import annotations
from typing import Literal

import mpmath as mp

from .high_precision_refinement import (
    mp_all_polynomial_roots,
    mp_p14_coefficients_from_ek,
)
from .parameter_sweep import lambda_from_ek
from .polynomial import candidate_roots, effective_polynomial_from_parameters

PolynomialBackend = Literal["float64", "mp"]


def polynomial_roots_from_ek(
    Arho: float | str,
    Amu: float | str,
    AG: float | str,
    Ek: float | str,
    *,
    backend: PolynomialBackend = "float64",
    decimal_precision: int = 80,
) -> tuple[list[complex] | list[mp.mpc], float | mp.mpf]:
    """Evaluate and solve P14 using the selected numerical backend.

    The default path is exactly the existing production implementation.  The
    ``mp`` path computes Lambda from Ek and evaluates every canonical
    coefficient without an intermediate binary64 Lambda.
    """

    if backend == "float64":
        lam = lambda_from_ek(float(Ek))
        effective = effective_polynomial_from_parameters(
            float(Arho), float(Amu), float(AG), lam
        )
        return [complex(root) for root in candidate_roots(effective)], lam
    if backend == "mp":
        coefficients, lam = mp_p14_coefficients_from_ek(
            Arho, Amu, AG, Ek, decimal_precision=decimal_precision
        )
        return mp_all_polynomial_roots(
            coefficients, decimal_precision=decimal_precision
        ), lam
    raise ValueError("backend must be 'float64' or 'mp'")
