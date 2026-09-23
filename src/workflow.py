"""Single-parameter-point orchestration for the finalized root workflow.

This module packages the production pipeline for one material-parameter tuple:

1. evaluate the verified polynomial coefficients;
2. compute polynomial candidate roots;
3. check the source-equation definition domain;
4. verify mathematically genuine roots against ``D_star_Lambda``;
5. apply the physical spatial-decay filter; and
6. assign one mutually exclusive parameter-point classification.

It is an orchestration layer only.  It does not duplicate coefficient
formulas, intermediate quantities, residual checks, or q-based admissibility
logic from the lower-level modules.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from .candidate_domain import CandidateDomainResult, check_candidates_domain
from .high_precision_refinement import (
    HighPrecisionOptions,
    refine_ambiguous_verification_run,
)
from .physical_admissibility import (
    PhysicalAdmissibilityRun,
    assess_physical_admissibility_many,
)
from .polynomial import (
    EffectivePolynomial,
    candidate_roots,
    effective_polynomial_from_parameters,
)
from .root_verification import VerificationRun, verify_candidate_roots


ParameterClassification = Literal[
    "outside_theorem_scope",
    "identically_zero_polynomial",
    "constant_nonzero_polynomial",
    "only_domain_invalid_candidates",
    "only_spurious_candidates",
    "genuine_roots_nonadmissible",
    "marginal_spatial_decay",
    "physically_admissible_roots",
    "numerically_unresolved",
]


TemporalStatus = Literal["growing", "decaying", "temporally_marginal"]


@dataclass(frozen=True)
class ParameterTuple:
    """Canonical nondimensional material-parameter tuple."""

    Arho: float
    Amu: float
    AG: float
    Lambda: float


@dataclass(frozen=True)
class TheoremScopeStatus:
    """Whether the parameter-level symbolic theorem scope applies.

    Candidate-level theorem conditions such as ``s != 0``, ``E_± != 0``, and
    ``B_± != 0`` are checked separately for each polynomial root in
    ``src.candidate_domain``.
    """

    applies: bool
    failed_conditions: tuple[str, ...]
    analysis_skipped: bool
    diagnostic_only: bool


@dataclass(frozen=True)
class WorkflowCounts:
    """Bookkeeping counts across the full one-point workflow."""

    effective_polynomial_degree: int | None
    total_polynomial_candidates: int
    domain_valid_candidates: int
    domain_invalid_candidates: int
    mathematically_genuine_roots: int
    ambiguous_mathematical_candidates: int
    physically_admissible_roots: int
    marginal_spatial_decay_roots: int
    physically_nonadmissible_roots: int
    evaluation_failures: int


@dataclass(frozen=True)
class PhysicalRootMetadata:
    """Temporal metadata for a physically admissible root.

    This is descriptive metadata only.  It is not used for accepting/rejecting
    roots and is not a physical-admissibility condition.
    """

    root: complex
    temporal_status: TemporalStatus
    real_part: float
    imaginary_part: float


@dataclass(frozen=True)
class ParameterPointResult:
    """Structured result for one parameter-point analysis."""

    parameters: ParameterTuple
    theorem_scope: TheoremScopeStatus
    effective_polynomial: EffectivePolynomial | None
    polynomial_candidates: tuple[complex, ...]
    domain_results: tuple[CandidateDomainResult, ...]
    verification_run: VerificationRun | None
    physical_run: PhysicalAdmissibilityRun | None
    counts: WorkflowCounts
    classification: ParameterClassification
    reduced_degree: bool
    warnings: tuple[str, ...]
    physical_root_metadata: tuple[PhysicalRootMetadata, ...]


def theorem_scope_status(
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    diagnostic_outside_scope: bool = False,
) -> TheoremScopeStatus:
    """Return the currently proven parameter-level theorem-scope status.

    The parameter-level part of the no-root-loss theorem is established for

    - finite real inputs;
    - ``Lambda > 0``;
    - ``-1 < Arho < 1``;
    - ``-1 < Amu < 1``;
    - ``-1 <= AG <= 1``.

    The limiting case ``Lambda = 0`` is outside the current parameter-level
    theorem scope and must not be silently treated as polynomial-complete.
    Candidate-level conditions are handled later by ``candidate_domain.py``.
    """

    failed: list[str] = []
    values = {
        "Arho": Arho,
        "Amu": Amu,
        "AG": AG,
        "Lambda": Lambda,
    }
    for name, value in values.items():
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            failed.append(f"{name}_finite_real")
            continue
        if not np.isfinite(numeric):
            failed.append(f"{name}_finite_real")

    if not failed:
        if not (-1.0 < float(Arho) < 1.0):
            failed.append("Arho_in_open_interval")
        if not (-1.0 < float(Amu) < 1.0):
            failed.append("Amu_in_open_interval")
        if not (-1.0 <= float(AG) <= 1.0):
            failed.append("AG_in_closed_interval")
        if not (float(Lambda) > 0.0):
            failed.append("Lambda_positive")

    applies = len(failed) == 0
    return TheoremScopeStatus(
        applies=applies,
        failed_conditions=tuple(failed),
        analysis_skipped=(not applies and not diagnostic_outside_scope),
        diagnostic_only=(not applies and diagnostic_outside_scope),
    )


def _empty_counts(effective_degree: int | None = None) -> WorkflowCounts:
    return WorkflowCounts(
        effective_polynomial_degree=effective_degree,
        total_polynomial_candidates=0,
        domain_valid_candidates=0,
        domain_invalid_candidates=0,
        mathematically_genuine_roots=0,
        ambiguous_mathematical_candidates=0,
        physically_admissible_roots=0,
        marginal_spatial_decay_roots=0,
        physically_nonadmissible_roots=0,
        evaluation_failures=0,
    )


def _temporal_status(root: complex, tolerance: float) -> TemporalStatus:
    if root.real > tolerance:
        return "growing"
    if root.real < -tolerance:
        return "decaying"
    return "temporally_marginal"


def _physical_metadata(
    physical_run: PhysicalAdmissibilityRun | None,
    *,
    temporal_tolerance: float,
) -> tuple[PhysicalRootMetadata, ...]:
    if physical_run is None:
        return ()
    rows: list[PhysicalRootMetadata] = []
    for result in physical_run.results:
        if result.physically_admissible:
            rows.append(
                PhysicalRootMetadata(
                    root=result.candidate,
                    temporal_status=_temporal_status(
                        result.candidate,
                        temporal_tolerance,
                    ),
                    real_part=float(result.candidate.real),
                    imaginary_part=float(result.candidate.imag),
                )
            )
    return tuple(rows)


def _counts(
    effective: EffectivePolynomial | None,
    candidates: tuple[complex, ...],
    domain_results: tuple[CandidateDomainResult, ...],
    verification_run: VerificationRun | None,
    physical_run: PhysicalAdmissibilityRun | None,
) -> WorkflowCounts:
    if verification_run is None or physical_run is None:
        return _empty_counts(None if effective is None else effective.effective_degree)

    verification = verification_run.summary
    physical = physical_run.summary
    return WorkflowCounts(
        effective_polynomial_degree=effective.effective_degree if effective else None,
        total_polynomial_candidates=len(candidates),
        domain_valid_candidates=verification.domain_valid_candidates,
        domain_invalid_candidates=verification.domain_invalid_candidates,
        mathematically_genuine_roots=verification.mathematically_genuine_roots,
        ambiguous_mathematical_candidates=verification.ambiguous_candidates,
        physically_admissible_roots=physical.physically_admissible_roots,
        marginal_spatial_decay_roots=physical.marginal_roots,
        physically_nonadmissible_roots=physical.physically_nonadmissible_roots,
        evaluation_failures=verification.evaluation_failed_candidates,
    )


def classify_parameter_point(
    theorem_scope: TheoremScopeStatus,
    effective: EffectivePolynomial | None,
    counts: WorkflowCounts,
    physical_run: PhysicalAdmissibilityRun | None,
    warnings: tuple[str, ...] = (),
) -> ParameterClassification:
    """Assign one mutually exclusive final classification."""

    if not theorem_scope.applies:
        return "outside_theorem_scope"
    if effective is None:
        return "numerically_unresolved"
    if effective.status == "identically_zero":
        return "identically_zero_polynomial"
    if effective.status == "constant_nonzero":
        return "constant_nonzero_polynomial"
    if (
        counts.ambiguous_mathematical_candidates > 0
        or counts.evaluation_failures > 0
        or "candidate_count_below_effective_degree" in warnings
    ):
        return "numerically_unresolved"
    if physical_run is not None and physical_run.summary.missing_domain_intermediates > 0:
        return "numerically_unresolved"
    if counts.total_polynomial_candidates == 0:
        return "numerically_unresolved"
    if counts.domain_valid_candidates == 0:
        return "only_domain_invalid_candidates"
    if counts.mathematically_genuine_roots == 0:
        return "only_spurious_candidates"
    if counts.physically_admissible_roots > 0:
        return "physically_admissible_roots"
    if counts.marginal_spatial_decay_roots > 0:
        return "marginal_spatial_decay"
    return "genuine_roots_nonadmissible"


def analyze_parameter_point(
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    analyze_outside_theorem_scope: bool = False,
    polynomial_relative_tolerance: float = 1e-12,
    polynomial_absolute_tolerance: float = 0.0,
    domain_relative_tolerance: float = 1e-10,
    domain_absolute_tolerance: float = 1e-14,
    residual_relative_tolerance: float = 1e-8,
    residual_absolute_tolerance: float = 1e-10,
    ambiguity_factor: float = 100.0,
    enable_high_precision_refinement: bool = False,
    high_precision_initial_decimal_precision: int = 80,
    high_precision_maximum_decimal_precision: int = 120,
    high_precision_precision_step: int = 20,
    high_precision_maximum_refinement_attempts: int = 80,
    high_precision_cluster_spacing_threshold: float = 1e-3,
    q_relative_tolerance: float = 1e-10,
    q_absolute_tolerance: float = 1e-12,
    temporal_tolerance: float = 1e-10,
) -> ParameterPointResult:
    """Run the complete single-parameter-point root workflow."""

    parameters = ParameterTuple(float(Arho), float(Amu), float(AG), float(Lambda))
    scope = theorem_scope_status(
        parameters.Arho,
        parameters.Amu,
        parameters.AG,
        parameters.Lambda,
        diagnostic_outside_scope=analyze_outside_theorem_scope,
    )
    warnings: list[str] = []
    if scope.diagnostic_only:
        warnings.append("analysis_outside_theorem_scope_is_diagnostic_only")

    if scope.analysis_skipped:
        return ParameterPointResult(
            parameters=parameters,
            theorem_scope=scope,
            effective_polynomial=None,
            polynomial_candidates=(),
            domain_results=(),
            verification_run=None,
            physical_run=None,
            counts=_empty_counts(),
            classification="outside_theorem_scope",
            reduced_degree=False,
            warnings=tuple(warnings),
            physical_root_metadata=(),
        )

    effective = effective_polynomial_from_parameters(
        parameters.Arho,
        parameters.Amu,
        parameters.AG,
        parameters.Lambda,
        relative_tolerance=polynomial_relative_tolerance,
        absolute_tolerance=polynomial_absolute_tolerance,
    )
    reduced_degree = effective.status == "reduced_degree"

    if effective.status in {"identically_zero", "constant_nonzero"}:
        counts = _empty_counts(effective.effective_degree)
        classification = classify_parameter_point(scope, effective, counts, None, tuple(warnings))
        return ParameterPointResult(
            parameters=parameters,
            theorem_scope=scope,
            effective_polynomial=effective,
            polynomial_candidates=(),
            domain_results=(),
            verification_run=None,
            physical_run=None,
            counts=counts,
            classification=classification,
            reduced_degree=reduced_degree,
            warnings=tuple(warnings),
            physical_root_metadata=(),
        )

    roots_array = candidate_roots(effective)
    candidates = tuple(complex(root) for root in roots_array)
    if effective.effective_degree is not None and len(candidates) < effective.effective_degree:
        warnings.append("candidate_count_below_effective_degree")

    domain_results = tuple(
        check_candidates_domain(
            candidates,
            parameters.Arho,
            parameters.Amu,
            parameters.AG,
            parameters.Lambda,
            relative_tolerance=domain_relative_tolerance,
            absolute_tolerance=domain_absolute_tolerance,
        )
    )
    verification_run = verify_candidate_roots(
        candidates,
        parameters.Arho,
        parameters.Amu,
        parameters.AG,
        parameters.Lambda,
        domain_results=domain_results,
        domain_relative_tolerance=domain_relative_tolerance,
        domain_absolute_tolerance=domain_absolute_tolerance,
        residual_relative_tolerance=residual_relative_tolerance,
        residual_absolute_tolerance=residual_absolute_tolerance,
        ambiguity_factor=ambiguity_factor,
    )
    if (
        enable_high_precision_refinement
        and verification_run.summary.ambiguous_candidates > 0
    ):
        verification_run = refine_ambiguous_verification_run(
            verification_run,
            effective,
            parameters.Arho,
            parameters.Amu,
            parameters.AG,
            parameters.Lambda,
            residual_relative_tolerance=residual_relative_tolerance,
            residual_absolute_tolerance=residual_absolute_tolerance,
            ambiguity_factor=ambiguity_factor,
            options=HighPrecisionOptions(
                initial_decimal_precision=high_precision_initial_decimal_precision,
                maximum_decimal_precision=high_precision_maximum_decimal_precision,
                precision_step=high_precision_precision_step,
                maximum_refinement_attempts=high_precision_maximum_refinement_attempts,
                cluster_spacing_threshold=high_precision_cluster_spacing_threshold,
            ),
            domain_relative_tolerance=domain_relative_tolerance,
            domain_absolute_tolerance=domain_absolute_tolerance,
        )
    physical_run = assess_physical_admissibility_many(
        verification_run.results,
        relative_tolerance=q_relative_tolerance,
        absolute_tolerance=q_absolute_tolerance,
    )
    counts = _counts(effective, candidates, domain_results, verification_run, physical_run)
    classification = classify_parameter_point(
        scope,
        effective,
        counts,
        physical_run,
        tuple(warnings),
    )

    return ParameterPointResult(
        parameters=parameters,
        theorem_scope=scope,
        effective_polynomial=effective,
        polynomial_candidates=candidates,
        domain_results=domain_results,
        verification_run=verification_run,
        physical_run=physical_run,
        counts=counts,
        classification=classification,
        reduced_degree=reduced_degree,
        warnings=tuple(warnings),
        physical_root_metadata=_physical_metadata(
            physical_run,
            temporal_tolerance=temporal_tolerance,
        ),
    )
