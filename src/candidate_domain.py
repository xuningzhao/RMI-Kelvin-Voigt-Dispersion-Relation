"""Definition-domain checks for P14 polynomial candidates.

This module evaluates the intermediate quantities that define the original
nondimensional dispersion relation on the principal square-root branch.  It
does not evaluate the dispersion residual and does not apply the physical
spatial-decay filter.

Notation follows ``THEORY.md#root-verification``:

    E_+ = s(1 + A_mu) + Lambda(1 + A_G)
    E_- = s(1 - A_mu) + Lambda(1 - A_G)

    R_+ = 1 + (1 + A_rho) s^2 / E_+
    R_- = 1 + (1 - A_rho) s^2 / E_-

    q_+ = sqrt(R_+),    q_- = sqrt(R_-)

    B_- = E_+ + E_- q_-
    B_+ = E_- + E_+ q_+

The square root is NumPy's principal complex square root, matching
``src.dispersion``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
from numpy.typing import ArrayLike


DEFAULT_RELATIVE_TOLERANCE = 1e-10
DEFAULT_ABSOLUTE_TOLERANCE = 1e-14


@dataclass(frozen=True)
class CandidateIntermediates:
    """Intermediate quantities used by the source dispersion relation."""

    s: complex
    E_plus: complex
    E_minus: complex
    R_plus: complex
    R_minus: complex
    q_plus: complex
    q_minus: complex
    B_plus: complex
    B_minus: complex


@dataclass(frozen=True)
class NonzeroCheck:
    """One scale-aware numerical nonzero check."""

    name: str
    value: complex
    magnitude: float
    threshold: float
    passed: bool


@dataclass(frozen=True)
class CandidateDomainResult:
    """Definition-domain diagnostic for one polynomial candidate."""

    s: complex
    domain_valid: bool
    failed_conditions: tuple[str, ...]
    intermediates: CandidateIntermediates | None
    checks: tuple[NonzeroCheck, ...]
    notes: tuple[str, ...] = ()


def _finite_complex(value: complex, name: str) -> complex:
    result = complex(value)
    if not np.isfinite(result.real) or not np.isfinite(result.imag):
        raise ValueError(f"{name} must be finite")
    return result


def _finite_real(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _threshold(
    *scale_terms: complex | float,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> float:
    scale = max([1.0, *[abs(complex(term)) for term in scale_terms]])
    return float(absolute_tolerance + relative_tolerance * scale)


def _nonzero_check(
    name: str,
    value: complex,
    *scale_terms: complex | float,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> NonzeroCheck:
    magnitude = float(abs(value))
    threshold = _threshold(
        *scale_terms,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )
    return NonzeroCheck(
        name=name,
        value=complex(value),
        magnitude=magnitude,
        threshold=threshold,
        passed=magnitude > threshold,
    )


def candidate_intermediates(
    s: complex,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
) -> CandidateIntermediates:
    """Compute ``E_±``, ``R_±``, principal ``q_±``, and ``B_±``.

    The caller is responsible for checking whether denominators are safely
    nonzero before interpreting quantities downstream of them.
    """

    s = _finite_complex(s, "s")
    Arho = _finite_real(Arho, "Arho")
    Amu = _finite_real(Amu, "Amu")
    AG = _finite_real(AG, "AG")
    Lambda = _finite_real(Lambda, "Lambda")

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        E_plus = s * (1.0 + Amu) + Lambda * (1.0 + AG)
        E_minus = s * (1.0 - Amu) + Lambda * (1.0 - AG)
        R_plus = 1.0 + ((1.0 + Arho) * s**2) / E_plus
        R_minus = 1.0 + ((1.0 - Arho) * s**2) / E_minus
        q_plus = np.sqrt(np.asarray(R_plus, dtype=np.complex128)).item()
        q_minus = np.sqrt(np.asarray(R_minus, dtype=np.complex128)).item()
        B_minus = E_plus + E_minus * q_minus
        B_plus = E_minus + E_plus * q_plus

    return CandidateIntermediates(
        s=complex(s),
        E_plus=complex(E_plus),
        E_minus=complex(E_minus),
        R_plus=complex(R_plus),
        R_minus=complex(R_minus),
        q_plus=complex(q_plus),
        q_minus=complex(q_minus),
        B_plus=complex(B_plus),
        B_minus=complex(B_minus),
    )


def check_candidate_domain(
    s: complex,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    relative_tolerance: float = DEFAULT_RELATIVE_TOLERANCE,
    absolute_tolerance: float = DEFAULT_ABSOLUTE_TOLERANCE,
) -> CandidateDomainResult:
    """Check the recorded nonzero definition conditions for one candidate."""

    if relative_tolerance < 0 or absolute_tolerance < 0:
        raise ValueError("domain tolerances must be nonnegative")

    s = _finite_complex(s, "s")
    Arho = _finite_real(Arho, "Arho")
    Amu = _finite_real(Amu, "Amu")
    AG = _finite_real(AG, "AG")
    Lambda = _finite_real(Lambda, "Lambda")

    checks: list[NonzeroCheck] = [
        _nonzero_check(
            "s",
            s,
            s,
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )
    ]

    notes: list[str] = []
    intermediates: CandidateIntermediates | None = None

    E_plus = s * (1.0 + Amu) + Lambda * (1.0 + AG)
    E_minus = s * (1.0 - Amu) + Lambda * (1.0 - AG)

    checks.append(
        _nonzero_check(
            "E_plus",
            E_plus,
            s * (1.0 + Amu),
            Lambda * (1.0 + AG),
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )
    )
    checks.append(
        _nonzero_check(
            "E_minus",
            E_minus,
            s * (1.0 - Amu),
            Lambda * (1.0 - AG),
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )
    )

    if checks[1].passed and checks[2].passed:
        intermediates = candidate_intermediates(s, Arho, Amu, AG, Lambda)
        if not all(
            np.isfinite(value.real) and np.isfinite(value.imag)
            for value in (
                intermediates.R_plus,
                intermediates.R_minus,
                intermediates.q_plus,
                intermediates.q_minus,
                intermediates.B_plus,
                intermediates.B_minus,
            )
        ):
            notes.append("nonfinite_intermediate")
        else:
            checks.append(
                _nonzero_check(
                    "B_plus",
                    intermediates.B_plus,
                    intermediates.E_minus,
                    intermediates.E_plus * intermediates.q_plus,
                    relative_tolerance=relative_tolerance,
                    absolute_tolerance=absolute_tolerance,
                )
            )
            checks.append(
                _nonzero_check(
                    "B_minus",
                    intermediates.B_minus,
                    intermediates.E_plus,
                    intermediates.E_minus * intermediates.q_minus,
                    relative_tolerance=relative_tolerance,
                    absolute_tolerance=absolute_tolerance,
                )
            )
    else:
        notes.append("skipped_radicals_after_E_denominator_failure")

    failed = tuple(check.name for check in checks if not check.passed)
    if "nonfinite_intermediate" in notes:
        failed = (*failed, "finite_intermediates")

    return CandidateDomainResult(
        s=s,
        domain_valid=len(failed) == 0,
        failed_conditions=failed,
        intermediates=intermediates,
        checks=tuple(checks),
        notes=tuple(notes),
    )


def check_candidates_domain(
    roots: Iterable[complex] | ArrayLike,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    relative_tolerance: float = DEFAULT_RELATIVE_TOLERANCE,
    absolute_tolerance: float = DEFAULT_ABSOLUTE_TOLERANCE,
) -> list[CandidateDomainResult]:
    """Check definition-domain conditions for many candidate roots."""

    return [
        check_candidate_domain(
            complex(root),
            Arho,
            Amu,
            AG,
            Lambda,
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )
        for root in np.asarray(list(roots), dtype=np.complex128).ravel()
    ]
