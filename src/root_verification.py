"""Mathematical verification of polynomial candidate roots.

Sprint 3 responsibility only:

1. accept polynomial candidates;
2. reuse definition-domain diagnostics;
3. evaluate the source nondimensional dispersion relation for domain-valid
   candidates; and
4. classify candidates as mathematically genuine or rejected/ambiguous.

This module determines *mathematical genuineness only*.  A mathematically
genuine root is a domain-valid polynomial candidate whose source-equation
residual passes the documented numerical tolerance.  Physical admissibility is
handled separately by the later spatial-decay filter based on ``q_+`` and
``q_-``.  This module deliberately does not apply physical admissibility
conditions, does not inspect ``Re(s)``, and does not perform parameter sweeps.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Literal

import numpy as np
from numpy.typing import ArrayLike

from .candidate_domain import (
    CandidateDomainResult,
    check_candidate_domain,
)
from .dispersion import D_star_Lambda
from .polynomial import candidate_roots, effective_polynomial_from_parameters


VerificationStatus = Literal[
    "mathematically_genuine",
    "domain_invalid",
    "residual_too_large",
    "numerically_ambiguous",
    "evaluation_failed",
]


DEFAULT_RESIDUAL_RELATIVE_TOLERANCE = 1e-8
DEFAULT_RESIDUAL_ABSOLUTE_TOLERANCE = 1e-10
DEFAULT_AMBIGUITY_FACTOR = 100.0


@dataclass(frozen=True)
class CandidateVerificationResult:
    """Mathematical verification record for one polynomial candidate."""

    candidate: complex
    domain_result: CandidateDomainResult
    domain_valid: bool
    absolute_residual: float | None
    relative_residual: float | None
    residual_scale: float | None
    verification_threshold: float | None
    mathematically_genuine: bool
    status: VerificationStatus
    rejection_reason: str | None
    dispersion_value: complex | None = None
    high_precision_refinement: Any | None = None


@dataclass(frozen=True)
class VerificationSummary:
    """Counts summarizing mathematical verification for one parameter tuple."""

    total_polynomial_candidates: int
    domain_valid_candidates: int
    mathematically_genuine_roots: int
    rejected_candidates: int
    ambiguous_candidates: int
    domain_invalid_candidates: int
    evaluation_failed_candidates: int


@dataclass(frozen=True)
class VerificationRun:
    """Full mathematical verification output for one parameter tuple."""

    results: tuple[CandidateVerificationResult, ...]
    summary: VerificationSummary

    @property
    def genuine_roots(self) -> tuple[complex, ...]:
        """Return candidates accepted as mathematically genuine roots."""

        return tuple(
            result.candidate
            for result in self.results
            if result.mathematically_genuine
        )


def residual_scale_from_domain(domain_result: CandidateDomainResult) -> float:
    """Return a scale for the source equation using stored intermediates.

    For domain-valid points, the source relation can be evaluated as

        D* = s**2 / B_- + s**2 / B_+ + 2.

    The residual scale is therefore the sum of magnitudes of those terms.
    A floor of one is applied later for relative residual reporting.
    """

    if not domain_result.domain_valid or domain_result.intermediates is None:
        raise ValueError("residual scale requires a domain-valid candidate")

    data = domain_result.intermediates
    with np.errstate(divide="raise", invalid="raise", over="raise"):
        term_minus = data.s**2 / data.B_minus
        term_plus = data.s**2 / data.B_plus
    scale = abs(term_minus) + abs(term_plus) + 2.0
    if not np.isfinite(scale):
        raise FloatingPointError("non-finite residual scale")
    return float(scale)


def verification_threshold(
    residual_scale: float,
    *,
    relative_tolerance: float = DEFAULT_RESIDUAL_RELATIVE_TOLERANCE,
    absolute_tolerance: float = DEFAULT_RESIDUAL_ABSOLUTE_TOLERANCE,
) -> float:
    """Return the absolute residual threshold for a scaled equation."""

    if residual_scale < 0 or not np.isfinite(residual_scale):
        raise ValueError("residual_scale must be finite and nonnegative")
    if relative_tolerance < 0 or absolute_tolerance < 0:
        raise ValueError("residual tolerances must be nonnegative")
    return float(absolute_tolerance + relative_tolerance * max(1.0, residual_scale))


def verify_candidate_root(
    candidate: complex,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    domain_result: CandidateDomainResult | None = None,
    domain_relative_tolerance: float = 1e-10,
    domain_absolute_tolerance: float = 1e-14,
    residual_relative_tolerance: float = DEFAULT_RESIDUAL_RELATIVE_TOLERANCE,
    residual_absolute_tolerance: float = DEFAULT_RESIDUAL_ABSOLUTE_TOLERANCE,
    ambiguity_factor: float = DEFAULT_AMBIGUITY_FACTOR,
) -> CandidateVerificationResult:
    """Verify one polynomial candidate against the source dispersion relation."""

    if ambiguity_factor < 1.0:
        raise ValueError("ambiguity_factor must be at least one")

    candidate = complex(candidate)
    if domain_result is None:
        domain_result = check_candidate_domain(
            candidate,
            Arho,
            Amu,
            AG,
            Lambda,
            relative_tolerance=domain_relative_tolerance,
            absolute_tolerance=domain_absolute_tolerance,
        )

    if not domain_result.domain_valid:
        return CandidateVerificationResult(
            candidate=candidate,
            domain_result=domain_result,
            domain_valid=False,
            absolute_residual=None,
            relative_residual=None,
            residual_scale=None,
            verification_threshold=None,
            mathematically_genuine=False,
            status="domain_invalid",
            rejection_reason="domain failure: "
            + ",".join(domain_result.failed_conditions),
        )

    try:
        value = complex(D_star_Lambda(candidate, Arho, Amu, AG, Lambda))
        absolute_residual = float(abs(value))
        scale = residual_scale_from_domain(domain_result)
        relative_residual = float(absolute_residual / max(1.0, scale))
        threshold = verification_threshold(
            scale,
            relative_tolerance=residual_relative_tolerance,
            absolute_tolerance=residual_absolute_tolerance,
        )
    except Exception as exc:  # noqa: BLE001 - preserve diagnostics from evaluator
        return CandidateVerificationResult(
            candidate=candidate,
            domain_result=domain_result,
            domain_valid=True,
            absolute_residual=None,
            relative_residual=None,
            residual_scale=None,
            verification_threshold=None,
            mathematically_genuine=False,
            status="evaluation_failed",
            rejection_reason=f"dispersion evaluation failed: {exc}",
        )

    if absolute_residual <= threshold:
        return CandidateVerificationResult(
            candidate=candidate,
            domain_result=domain_result,
            domain_valid=True,
            absolute_residual=absolute_residual,
            relative_residual=relative_residual,
            residual_scale=scale,
            verification_threshold=threshold,
            mathematically_genuine=True,
            status="mathematically_genuine",
            rejection_reason=None,
            dispersion_value=value,
        )

    if absolute_residual <= ambiguity_factor * threshold:
        status: VerificationStatus = "numerically_ambiguous"
        reason = "residual near verification threshold"
    else:
        status = "residual_too_large"
        reason = "residual too large"

    return CandidateVerificationResult(
        candidate=candidate,
        domain_result=domain_result,
        domain_valid=True,
        absolute_residual=absolute_residual,
        relative_residual=relative_residual,
        residual_scale=scale,
        verification_threshold=threshold,
        mathematically_genuine=False,
        status=status,
        rejection_reason=reason,
        dispersion_value=value,
    )


def summarize_verification(
    results: Iterable[CandidateVerificationResult],
) -> VerificationSummary:
    """Summarize candidate verification records."""

    result_tuple = tuple(results)
    total = len(result_tuple)
    domain_valid = sum(result.domain_valid for result in result_tuple)
    genuine = sum(result.mathematically_genuine for result in result_tuple)
    ambiguous = sum(result.status == "numerically_ambiguous" for result in result_tuple)
    domain_invalid = sum(result.status == "domain_invalid" for result in result_tuple)
    evaluation_failed = sum(result.status == "evaluation_failed" for result in result_tuple)
    rejected = total - genuine

    return VerificationSummary(
        total_polynomial_candidates=total,
        domain_valid_candidates=domain_valid,
        mathematically_genuine_roots=genuine,
        rejected_candidates=rejected,
        ambiguous_candidates=ambiguous,
        domain_invalid_candidates=domain_invalid,
        evaluation_failed_candidates=evaluation_failed,
    )


def verify_candidate_roots(
    candidates: Iterable[complex] | ArrayLike,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    domain_results: Iterable[CandidateDomainResult] | None = None,
    domain_relative_tolerance: float = 1e-10,
    domain_absolute_tolerance: float = 1e-14,
    residual_relative_tolerance: float = DEFAULT_RESIDUAL_RELATIVE_TOLERANCE,
    residual_absolute_tolerance: float = DEFAULT_RESIDUAL_ABSOLUTE_TOLERANCE,
    ambiguity_factor: float = DEFAULT_AMBIGUITY_FACTOR,
) -> VerificationRun:
    """Verify many polynomial candidates against ``D_star_Lambda``."""

    roots = [complex(root) for root in np.asarray(list(candidates), dtype=np.complex128).ravel()]
    domain_list = None if domain_results is None else tuple(domain_results)
    if domain_list is not None and len(domain_list) != len(roots):
        raise ValueError("domain_results length must match candidates length")

    results = []
    for index, root in enumerate(roots):
        results.append(
            verify_candidate_root(
                root,
                Arho,
                Amu,
                AG,
                Lambda,
                domain_result=None if domain_list is None else domain_list[index],
                domain_relative_tolerance=domain_relative_tolerance,
                domain_absolute_tolerance=domain_absolute_tolerance,
                residual_relative_tolerance=residual_relative_tolerance,
                residual_absolute_tolerance=residual_absolute_tolerance,
                ambiguity_factor=ambiguity_factor,
            )
        )

    result_tuple = tuple(results)
    return VerificationRun(
        results=result_tuple,
        summary=summarize_verification(result_tuple),
    )


def verify_p14_roots_from_parameters(
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    polynomial_relative_tolerance: float = 1e-12,
    polynomial_absolute_tolerance: float = 0.0,
    domain_relative_tolerance: float = 1e-10,
    domain_absolute_tolerance: float = 1e-14,
    residual_relative_tolerance: float = DEFAULT_RESIDUAL_RELATIVE_TOLERANCE,
    residual_absolute_tolerance: float = DEFAULT_RESIDUAL_ABSOLUTE_TOLERANCE,
    ambiguity_factor: float = DEFAULT_AMBIGUITY_FACTOR,
) -> VerificationRun:
    """Generate P14 candidates and verify the mathematically genuine subset."""

    effective = effective_polynomial_from_parameters(
        Arho,
        Amu,
        AG,
        Lambda,
        relative_tolerance=polynomial_relative_tolerance,
        absolute_tolerance=polynomial_absolute_tolerance,
    )
    roots = candidate_roots(effective)
    return verify_candidate_roots(
        roots,
        Arho,
        Amu,
        AG,
        Lambda,
        domain_relative_tolerance=domain_relative_tolerance,
        domain_absolute_tolerance=domain_absolute_tolerance,
        residual_relative_tolerance=residual_relative_tolerance,
        residual_absolute_tolerance=residual_absolute_tolerance,
        ambiguity_factor=ambiguity_factor,
    )
