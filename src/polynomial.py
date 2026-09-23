"""Polynomial candidate infrastructure for the verified P14 equation.

This module handles only Sprint 1 responsibilities:

- evaluate the verified coefficient list;
- detect the structurally effective degree;
- normalize coefficients for numerical root computation;
- compute finite complex roots of the effective polynomial; and
- optionally cluster numerically duplicate roots while retaining multiplicity.

It deliberately does *not* verify candidates against the original dispersion
relation and does *not* apply physical admissibility filters.  Those belong to
later workflow stages.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .polynomial_coefficients import GENERIC_DEGREE, p14_coefficients


PolynomialStatus = Literal[
    "generic_degree",
    "reduced_degree",
    "near_structural_degeneracy",
    "constant_nonzero",
    "identically_zero",
]


STRUCTURAL_ABSOLUTE_TOLERANCE = 10.0 * np.finfo(float).eps
"""Absolute floor for recognizing floating-point structural zeros.

Degree reduction for the verified P14 polynomial must be theorem-aware, not
based on the largest coefficient in a badly scaled polynomial.  This tolerance
is used only to absorb floating-point representation noise in coefficients
that are structurally zero, for example when ``A_rho + A_mu == 0``.
"""


@dataclass(frozen=True)
class EffectivePolynomial:
    """Numerically effective representation of ``P14`` at one parameter tuple.

    Attributes
    ----------
    raw_coefficients:
        Low-to-high coefficients ``[c0, ..., c14]``.
    coefficient_scale:
        Maximum absolute raw coefficient.  Used for diagnostics and scalar
        normalization, not for deleting leading coefficients.
    trim_tolerance:
        Absolute structural-zero threshold applied to coefficient magnitudes.
    effective_degree:
        Highest retained power of ``s``.  ``None`` for the identically-zero
        structural case.
    coefficients_ascending:
        Low-to-high retained coefficients.
    coefficients_descending:
        High-to-low retained coefficients, as expected by ``numpy.roots``.
    normalized_descending:
        High-to-low coefficients of the scaled-variable polynomial divided by
        ``normalization_factor``.  If ``s = alpha z``, these are coefficients
        for the polynomial in ``z``.
    normalization_factor:
        Maximum absolute retained coefficient.  Scalar normalization does not
        change polynomial roots but improves numerical scaling.
    root_scale:
        Positive scale ``alpha`` used in ``s = alpha z``.
    status:
        Structural degree category.
    """

    raw_coefficients: NDArray[np.complex128]
    coefficient_scale: float
    trim_tolerance: float
    effective_degree: int | None
    coefficients_ascending: NDArray[np.complex128]
    coefficients_descending: NDArray[np.complex128]
    normalized_descending: NDArray[np.complex128]
    normalization_factor: float
    root_scale: float
    status: PolynomialStatus


@dataclass(frozen=True)
class RootCluster:
    """Cluster of numerically nearby candidate roots.

    ``multiplicity`` is the number of raw numerical roots in the cluster; it is
    not a symbolic algebraic multiplicity proof.
    """

    root: complex
    multiplicity: int
    members: tuple[complex, ...]
    radius: float


def evaluate_coefficients(
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
) -> NDArray[np.complex128]:
    """Evaluate verified low-to-high coefficients ``[c0, ..., c14]``."""

    return np.asarray(
        p14_coefficients(Arho, Amu, AG, Lambda),
        dtype=np.complex128,
    )


def polynomial_value(s: complex, coefficients_ascending: ArrayLike) -> complex:
    """Evaluate a low-to-high coefficient polynomial at complex ``s``."""

    coeffs = np.asarray(coefficients_ascending, dtype=np.complex128)
    z = np.asarray(s, dtype=np.complex128)
    value = np.asarray(0.0, dtype=np.complex128)
    for coeff in coeffs[::-1]:
        value = value * z + coeff
    return complex(value.item())


def derivative_coefficients(coefficients_ascending: ArrayLike) -> NDArray[np.complex128]:
    """Return low-to-high coefficients of the polynomial derivative."""

    coeffs = np.asarray(coefficients_ascending, dtype=np.complex128)
    if coeffs.size <= 1:
        return np.asarray([0.0], dtype=np.complex128)
    powers = np.arange(1, coeffs.size, dtype=np.float64)
    return coeffs[1:] * powers


def _structural_tolerance(absolute_tolerance: float) -> float:
    return float(max(absolute_tolerance, STRUCTURAL_ABSOLUTE_TOLERANCE))


def _scaled_variable_coefficients(
    coefficients_ascending: NDArray[np.complex128],
    root_scale: float,
) -> NDArray[np.complex128]:
    powers = np.arange(coefficients_ascending.size, dtype=np.float64)
    return coefficients_ascending * np.power(root_scale, powers)


def choose_root_scale(coefficients_ascending: ArrayLike) -> float:
    """Choose a positive scale ``alpha`` for solving in ``s = alpha z``.

    The scale is the geometric mean of coefficient-ratio root-radius estimates

    ``(|c_j| / |c_n|)**(1 / (n - j))``

    over nonzero lower coefficients.  This balances the polynomial for the
    root solver without discarding mathematically essential leading terms.
    """

    coeffs = np.asarray(coefficients_ascending, dtype=np.complex128)
    if coeffs.ndim != 1:
        raise ValueError("coefficients must be a one-dimensional array")
    degree = coeffs.size - 1
    if degree <= 0:
        return 1.0
    leading = float(abs(coeffs[-1]))
    if leading == 0.0 or not np.isfinite(leading):
        return 1.0

    estimates: list[float] = []
    for power, coeff in enumerate(coeffs[:-1]):
        magnitude = float(abs(coeff))
        if magnitude == 0.0 or not np.isfinite(magnitude):
            continue
        estimate = (magnitude / leading) ** (1.0 / (degree - power))
        if np.isfinite(estimate) and estimate > 0.0:
            estimates.append(float(estimate))

    if not estimates:
        return 1.0
    log_scale = float(np.mean(np.log(estimates)))
    root_scale = float(np.exp(log_scale))
    if not np.isfinite(root_scale) or root_scale <= 0.0:
        return 1.0
    return root_scale


def polish_polynomial_roots(
    roots: ArrayLike,
    coefficients_ascending: ArrayLike,
    *,
    max_iterations: int = 20,
    step_tolerance: float = 1e-13,
    derivative_tolerance: float = 1e-14,
) -> NDArray[np.complex128]:
    """Polish candidate roots by Newton iteration on the polynomial itself.

    The method is intentionally conservative: if the derivative is too small
    or an iteration becomes non-finite, the last finite iterate is retained.
    This improves ordinary simple roots without pretending to solve symbolic
    multiplicity issues.
    """

    if max_iterations < 0:
        raise ValueError("max_iterations must be nonnegative")
    if step_tolerance < 0 or derivative_tolerance < 0:
        raise ValueError("polishing tolerances must be nonnegative")

    roots_array = np.asarray(roots, dtype=np.complex128).ravel()
    coeffs = np.asarray(coefficients_ascending, dtype=np.complex128)
    deriv = derivative_coefficients(coeffs)
    polished: list[complex] = []

    for root in roots_array:
        z = complex(root)
        for _ in range(max_iterations):
            value = polynomial_value(z, coeffs)
            slope = polynomial_value(z, deriv)
            if abs(slope) <= derivative_tolerance:
                break
            step = value / slope
            next_z = z - step
            if not np.isfinite(next_z.real) or not np.isfinite(next_z.imag):
                break
            z = next_z
            if abs(step) <= step_tolerance * max(1.0, abs(z)):
                break
        polished.append(z)

    return np.asarray(polished, dtype=np.complex128)


def build_effective_polynomial(
    coefficients_ascending: ArrayLike,
    *,
    relative_tolerance: float = 1e-12,
    absolute_tolerance: float = 0.0,
) -> EffectivePolynomial:
    """Build a structurally degree-aware normalized root-solver input.

    Degree reduction is based only on an absolute structural-zero tolerance,
    not on the global coefficient scale.  This is essential for P14: in
    high-``Lambda`` cases lower coefficients can be enormous while the true
    leading coefficient remains moderate.  Dropping that leading coefficient
    would lose roots before source-equation verification.

    The coefficient order of the returned ``coefficients_ascending`` remains
    low-to-high.  Only the root-solver arrays are descending.
    """

    if relative_tolerance < 0:
        raise ValueError("relative_tolerance must be nonnegative")
    if absolute_tolerance < 0:
        raise ValueError("absolute_tolerance must be nonnegative")

    raw = np.asarray(coefficients_ascending, dtype=np.complex128)
    if raw.ndim != 1:
        raise ValueError("coefficients must be a one-dimensional array")
    if raw.size == 0:
        raise ValueError("at least one coefficient is required")
    if not np.all(np.isfinite(raw.real) & np.isfinite(raw.imag)):
        raise ValueError("coefficients must be finite")

    coefficient_scale = float(np.max(np.abs(raw)))
    trim_tolerance = _structural_tolerance(absolute_tolerance)

    if coefficient_scale <= trim_tolerance:
        empty = np.asarray([], dtype=np.complex128)
        return EffectivePolynomial(
            raw_coefficients=raw,
            coefficient_scale=coefficient_scale,
            trim_tolerance=trim_tolerance,
            effective_degree=None,
            coefficients_ascending=empty,
            coefficients_descending=empty,
            normalized_descending=empty,
            normalization_factor=0.0,
            root_scale=1.0,
            status="identically_zero",
        )

    retained = np.flatnonzero(np.abs(raw) > trim_tolerance)
    effective_degree = int(retained[-1])
    trimmed_ascending = raw[: effective_degree + 1].copy()
    descending = trimmed_ascending[::-1].copy()
    root_scale = choose_root_scale(trimmed_ascending)
    scaled_ascending = _scaled_variable_coefficients(trimmed_ascending, root_scale)
    scaled_descending = scaled_ascending[::-1].copy()
    normalization_factor = float(np.max(np.abs(scaled_ascending)))
    normalized_descending = scaled_descending / normalization_factor

    if effective_degree == 0:
        status: PolynomialStatus = "constant_nonzero"
    elif effective_degree < GENERIC_DEGREE:
        status = "reduced_degree"
    elif abs(trimmed_ascending[-1]) <= 100.0 * trim_tolerance:
        status = "near_structural_degeneracy"
    else:
        status = "generic_degree"

    return EffectivePolynomial(
        raw_coefficients=raw,
        coefficient_scale=coefficient_scale,
        trim_tolerance=trim_tolerance,
        effective_degree=effective_degree,
        coefficients_ascending=trimmed_ascending,
        coefficients_descending=descending,
        normalized_descending=normalized_descending,
        normalization_factor=normalization_factor,
        root_scale=root_scale,
        status=status,
    )


def effective_polynomial_from_parameters(
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    relative_tolerance: float = 1e-12,
    absolute_tolerance: float = 0.0,
) -> EffectivePolynomial:
    """Evaluate structurally degree-aware ``P14`` coefficients for one tuple."""

    return build_effective_polynomial(
        evaluate_coefficients(Arho, Amu, AG, Lambda),
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )


def candidate_roots(
    effective: EffectivePolynomial,
    *,
    polish: bool = True,
    polish_max_iterations: int = 20,
) -> NDArray[np.complex128]:
    """Compute all finite complex roots of the effective polynomial."""

    if effective.effective_degree is None or effective.effective_degree <= 0:
        return np.asarray([], dtype=np.complex128)

    z_roots = np.roots(effective.normalized_descending)
    roots = np.asarray(z_roots, dtype=np.complex128) * effective.root_scale
    roots = np.asarray(roots, dtype=np.complex128)
    finite = np.isfinite(roots.real) & np.isfinite(roots.imag)
    roots = roots[finite]
    if polish:
        roots = polish_polynomial_roots(
            roots,
            effective.coefficients_ascending,
            max_iterations=polish_max_iterations,
        )
    return roots


def _near_root(a: complex, b: complex, absolute_tolerance: float, relative_tolerance: float) -> bool:
    scale = max(1.0, abs(a), abs(b))
    return abs(a - b) <= absolute_tolerance + relative_tolerance * scale


def merge_duplicate_roots(
    roots: ArrayLike,
    *,
    absolute_tolerance: float = 1e-8,
    relative_tolerance: float = 1e-8,
) -> list[RootCluster]:
    """Cluster nearby numerical roots while preserving multiplicity counts."""

    if absolute_tolerance < 0 or relative_tolerance < 0:
        raise ValueError("root-merge tolerances must be nonnegative")

    remaining = [complex(z) for z in np.asarray(roots, dtype=np.complex128).ravel()]
    clusters: list[RootCluster] = []

    while remaining:
        seed = remaining.pop(0)
        members = [seed]
        changed = True
        while changed:
            changed = False
            center = sum(members) / len(members)
            keep: list[complex] = []
            for root in remaining:
                if _near_root(root, center, absolute_tolerance, relative_tolerance):
                    members.append(root)
                    changed = True
                else:
                    keep.append(root)
            remaining = keep

        center = sum(members) / len(members)
        radius = max(abs(root - center) for root in members)
        clusters.append(
            RootCluster(
                root=complex(center),
                multiplicity=len(members),
                members=tuple(members),
                radius=float(radius),
            )
        )

    clusters.sort(key=lambda cluster: (cluster.root.real, cluster.root.imag))
    return clusters


def p14_candidate_roots(
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    relative_tolerance: float = 1e-12,
    absolute_tolerance: float = 0.0,
    merge_duplicates: bool = False,
    root_merge_absolute_tolerance: float = 1e-8,
    root_merge_relative_tolerance: float = 1e-8,
    polish_roots: bool = True,
) -> tuple[EffectivePolynomial, NDArray[np.complex128], list[RootCluster] | None]:
    """Convenience wrapper returning effective polynomial and candidate roots."""

    effective = effective_polynomial_from_parameters(
        Arho,
        Amu,
        AG,
        Lambda,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )
    roots = candidate_roots(effective, polish=polish_roots)
    clusters = (
        merge_duplicate_roots(
            roots,
            absolute_tolerance=root_merge_absolute_tolerance,
            relative_tolerance=root_merge_relative_tolerance,
        )
        if merge_duplicates
        else None
    )
    return effective, roots, clusters
