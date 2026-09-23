"""Physical admissibility filter for mathematically genuine roots.

Sprint 4 responsibility only:

1. accept mathematical verification records from ``src.root_verification``;
2. reuse the already-computed ``q_+`` and ``q_-`` values stored in the
   definition-domain diagnostics; and
3. classify the semi-infinite normal-mode spatial decay conditions.

This module does not compute polynomial candidates, does not evaluate
``D_star_Lambda``, does not inspect ``Re(s)``, and does not classify temporal
stability.  It only checks whether mathematically genuine roots satisfy

    Re(q_+) > 0,    Re(q_-) > 0.

With the principal complex square root used upstream, ``Re(q_±)`` is normally
nonnegative.  Therefore realistic boundary cases usually appear as marginal
spatial decay rather than strongly negative decay.  Negative real parts are
still handled explicitly as a robustness check against invalid or nonprincipal
upstream values.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np

from .root_verification import CandidateVerificationResult


PhysicalStatus = Literal[
    "physically_admissible",
    "marginal_spatial_decay",
    "physically_nonadmissible",
    "not_mathematically_genuine",
    "missing_domain_intermediates",
]


DEFAULT_Q_RELATIVE_TOLERANCE = 1e-10
DEFAULT_Q_ABSOLUTE_TOLERANCE = 1e-12


@dataclass(frozen=True)
class PhysicalAdmissibilityResult:
    """Physical spatial-decay diagnostic for one candidate."""

    candidate: complex
    mathematically_genuine: bool
    q_plus: complex | None
    q_minus: complex | None
    re_q_plus: float | None
    re_q_minus: float | None
    tolerance_plus: float | None
    tolerance_minus: float | None
    physically_admissible: bool
    status: PhysicalStatus
    reason: str
    verification_result: CandidateVerificationResult


@dataclass(frozen=True)
class PhysicalAdmissibilitySummary:
    """Counts summarizing physical admissibility for one candidate set."""

    total_candidates: int
    mathematically_genuine_roots: int
    physically_admissible_roots: int
    marginal_roots: int
    physically_nonadmissible_roots: int
    skipped_not_genuine: int
    missing_domain_intermediates: int


@dataclass(frozen=True)
class PhysicalAdmissibilityRun:
    """Full physical admissibility output for one candidate set."""

    results: tuple[PhysicalAdmissibilityResult, ...]
    summary: PhysicalAdmissibilitySummary

    @property
    def physically_admissible_roots(self) -> tuple[complex, ...]:
        """Return roots accepted by the spatial-decay filter."""

        return tuple(
            result.candidate
            for result in self.results
            if result.physically_admissible
        )


def q_real_tolerance(
    q: complex,
    *,
    relative_tolerance: float = DEFAULT_Q_RELATIVE_TOLERANCE,
    absolute_tolerance: float = DEFAULT_Q_ABSOLUTE_TOLERANCE,
) -> float:
    """Return a scale-aware tolerance for deciding ``Re(q) > 0``."""

    if relative_tolerance < 0 or absolute_tolerance < 0:
        raise ValueError("q-real tolerances must be nonnegative")
    q = complex(q)
    if not np.isfinite(q.real) or not np.isfinite(q.imag):
        raise ValueError("q must be finite")
    return float(absolute_tolerance + relative_tolerance * max(1.0, abs(q)))


def classify_q_decay(
    q_plus: complex,
    q_minus: complex,
    *,
    relative_tolerance: float = DEFAULT_Q_RELATIVE_TOLERANCE,
    absolute_tolerance: float = DEFAULT_Q_ABSOLUTE_TOLERANCE,
) -> tuple[PhysicalStatus, str, bool, float, float]:
    """Classify the spatial-decay inequalities for stored ``q_±`` values."""

    q_plus = complex(q_plus)
    q_minus = complex(q_minus)
    tol_plus = q_real_tolerance(
        q_plus,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )
    tol_minus = q_real_tolerance(
        q_minus,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )

    upper_positive = q_plus.real > tol_plus
    lower_positive = q_minus.real > tol_minus
    upper_marginal = abs(q_plus.real) <= tol_plus
    lower_marginal = abs(q_minus.real) <= tol_minus

    if upper_positive and lower_positive:
        return (
            "physically_admissible",
            "spatial_decay_satisfied",
            True,
            tol_plus,
            tol_minus,
        )

    if upper_marginal or lower_marginal:
        return (
            "marginal_spatial_decay",
            "marginal_spatial_decay",
            False,
            tol_plus,
            tol_minus,
        )

    failed_upper = not upper_positive
    failed_lower = not lower_positive
    if failed_upper and failed_lower:
        reason = "failed_both_decay"
    elif failed_upper:
        reason = "failed_upper_decay"
    else:
        reason = "failed_lower_decay"
    return (
        "physically_nonadmissible",
        reason,
        False,
        tol_plus,
        tol_minus,
    )


def assess_physical_admissibility(
    verification_result: CandidateVerificationResult,
    *,
    relative_tolerance: float = DEFAULT_Q_RELATIVE_TOLERANCE,
    absolute_tolerance: float = DEFAULT_Q_ABSOLUTE_TOLERANCE,
) -> PhysicalAdmissibilityResult:
    """Apply the spatial-decay filter to one verification record."""

    if not verification_result.mathematically_genuine:
        return PhysicalAdmissibilityResult(
            candidate=verification_result.candidate,
            mathematically_genuine=False,
            q_plus=None,
            q_minus=None,
            re_q_plus=None,
            re_q_minus=None,
            tolerance_plus=None,
            tolerance_minus=None,
            physically_admissible=False,
            status="not_mathematically_genuine",
            reason="not_mathematically_genuine",
            verification_result=verification_result,
        )

    intermediates = verification_result.domain_result.intermediates
    if intermediates is None:
        return PhysicalAdmissibilityResult(
            candidate=verification_result.candidate,
            mathematically_genuine=True,
            q_plus=None,
            q_minus=None,
            re_q_plus=None,
            re_q_minus=None,
            tolerance_plus=None,
            tolerance_minus=None,
            physically_admissible=False,
            status="missing_domain_intermediates",
            reason="missing_domain_intermediates",
            verification_result=verification_result,
        )

    status, reason, accepted, tol_plus, tol_minus = classify_q_decay(
        intermediates.q_plus,
        intermediates.q_minus,
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
    )
    return PhysicalAdmissibilityResult(
        candidate=verification_result.candidate,
        mathematically_genuine=True,
        q_plus=intermediates.q_plus,
        q_minus=intermediates.q_minus,
        re_q_plus=float(intermediates.q_plus.real),
        re_q_minus=float(intermediates.q_minus.real),
        tolerance_plus=tol_plus,
        tolerance_minus=tol_minus,
        physically_admissible=accepted,
        status=status,
        reason=reason,
        verification_result=verification_result,
    )


def summarize_physical_admissibility(
    results: Iterable[PhysicalAdmissibilityResult],
) -> PhysicalAdmissibilitySummary:
    """Summarize physical admissibility records."""

    result_tuple = tuple(results)
    return PhysicalAdmissibilitySummary(
        total_candidates=len(result_tuple),
        mathematically_genuine_roots=sum(
            result.mathematically_genuine for result in result_tuple
        ),
        physically_admissible_roots=sum(
            result.status == "physically_admissible" for result in result_tuple
        ),
        marginal_roots=sum(
            result.status == "marginal_spatial_decay" for result in result_tuple
        ),
        physically_nonadmissible_roots=sum(
            result.status == "physically_nonadmissible" for result in result_tuple
        ),
        skipped_not_genuine=sum(
            result.status == "not_mathematically_genuine" for result in result_tuple
        ),
        missing_domain_intermediates=sum(
            result.status == "missing_domain_intermediates" for result in result_tuple
        ),
    )


def assess_physical_admissibility_many(
    verification_results: Iterable[CandidateVerificationResult],
    *,
    relative_tolerance: float = DEFAULT_Q_RELATIVE_TOLERANCE,
    absolute_tolerance: float = DEFAULT_Q_ABSOLUTE_TOLERANCE,
) -> PhysicalAdmissibilityRun:
    """Apply physical admissibility checks to many verification records."""

    results = tuple(
        assess_physical_admissibility(
            result,
            relative_tolerance=relative_tolerance,
            absolute_tolerance=absolute_tolerance,
        )
        for result in verification_results
    )
    return PhysicalAdmissibilityRun(
        results=results,
        summary=summarize_physical_admissibility(results),
    )

