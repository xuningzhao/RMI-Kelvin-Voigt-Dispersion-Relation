"""Targeted high-precision refinement for ambiguous polynomial candidates.

This module is deliberately narrow.  The standard production workflow remains
double precision for ordinary candidates.  High precision is invoked only for
candidate roots whose mathematical verification status is
``numerically_ambiguous``.

The refinement recomputes the verified P14 coefficients with arbitrary
precision, polishes the nearby polynomial root, and evaluates the original
``D_star_Lambda`` expression with the same principal square-root convention.
It then resolves the candidate as mathematically genuine, spurious, or still
ambiguous without relaxing any global tolerances.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Iterable, Literal

import mpmath as mp
import numpy as np

from .candidate_domain import check_candidate_domain
from .polynomial import EffectivePolynomial
from .root_verification import (
    CandidateVerificationResult,
    VerificationRun,
    summarize_verification,
)


RefinementStatus = Literal[
    "not_attempted",
    "refined_genuine",
    "refined_spurious",
    "still_ambiguous",
    "refinement_failed",
]


@dataclass(frozen=True)
class HighPrecisionOptions:
    """Configuration for targeted high-precision refinement."""

    initial_decimal_precision: int = 80
    maximum_decimal_precision: int = 120
    precision_step: int = 20
    maximum_refinement_attempts: int = 80
    cluster_spacing_threshold: float = 1e-3
    polynomial_residual_tolerance: float = 1e-50


@dataclass(frozen=True)
class HighPrecisionRefinementResult:
    """Diagnostic record for one high-precision refinement attempt."""

    original_candidate: complex
    refined_candidate: complex | None
    working_precision: int
    high_precision_p14_residual: float | None
    high_precision_absolute_residual: float | None
    high_precision_relative_residual: float | None
    verification_threshold: float | None
    status: RefinementStatus
    reason: str
    cluster_size: int
    nearest_root_spacing: float | None
    full_polynomial_solve_used: bool


def _coefficient_list_ast() -> ast.List:
    source_path = Path(__file__).with_name("polynomial_coefficients.py")
    module = ast.parse(source_path.read_text(encoding="utf-8"))

    class Visitor(ast.NodeVisitor):
        def __init__(self) -> None:
            self.node: ast.List | None = None

        def visit_Call(self, node: ast.Call) -> None:  # noqa: N802 - ast API
            if (
                isinstance(node.func, ast.Attribute)
                and node.func.attr == "asarray"
                and node.args
                and isinstance(node.args[0], ast.List)
                and len(node.args[0].elts) == 15
            ):
                self.node = node.args[0]
            self.generic_visit(node)

    visitor = Visitor()
    visitor.visit(module)
    if visitor.node is None:
        raise RuntimeError("could not locate P14 coefficient list")
    return visitor.node


@lru_cache(maxsize=1)
def _compiled_coefficient_expressions() -> tuple[object, ...]:
    """Return compiled expressions from the canonical Python coefficient file."""

    expressions = []
    for node in _coefficient_list_ast().elts:
        expressions.append(compile(ast.Expression(node), "<P14 coefficient>", "eval"))
    return tuple(expressions)


def mp_p14_coefficients(
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    decimal_precision: int,
) -> list[mp.mpc]:
    """Evaluate low-to-high P14 coefficients at arbitrary precision."""

    mp.mp.dps = decimal_precision
    environment = {
        "Arho": mp.mpf(str(float(Arho))),
        "Amu": mp.mpf(str(float(Amu))),
        "AG": mp.mpf(str(float(AG))),
        "Lambda": mp.mpf(str(float(Lambda))),
    }
    return [
        mp.mpc(eval(expr, {"__builtins__": {}}, environment))
        for expr in _compiled_coefficient_expressions()
    ]


def mp_lambda_from_ek(Ek: str | float | mp.mpf, *, decimal_precision: int) -> mp.mpf:
    """Compute Lambda from Ek without an intermediate binary64 Lambda."""

    mp.mp.dps = decimal_precision
    ek = Ek if isinstance(Ek, mp.mpf) else mp.mpf(str(Ek))
    if not (-1 < ek < 1):
        raise ValueError("mp backend requires -1 < Ek < 1")
    return ((1 - ek) / (1 + ek)) ** 2


def mp_p14_coefficients_from_ek(
    Arho: str | float | mp.mpf,
    Amu: str | float | mp.mpf,
    AG: str | float | mp.mpf,
    Ek: str | float | mp.mpf,
    *,
    decimal_precision: int,
) -> tuple[list[mp.mpc], mp.mpf]:
    """Evaluate the canonical P14 expressions directly from an mp Ek value."""

    mp.mp.dps = decimal_precision
    environment = {
        "Arho": Arho if isinstance(Arho, mp.mpf) else mp.mpf(str(Arho)),
        "Amu": Amu if isinstance(Amu, mp.mpf) else mp.mpf(str(Amu)),
        "AG": AG if isinstance(AG, mp.mpf) else mp.mpf(str(AG)),
        "Lambda": mp_lambda_from_ek(Ek, decimal_precision=decimal_precision),
    }
    coefficients = [
        mp.mpc(eval(expr, {"__builtins__": {}}, environment))
        for expr in _compiled_coefficient_expressions()
    ]
    return coefficients, environment["Lambda"]


def mp_all_polynomial_roots(
    coefficients: list[mp.mpc],
    *,
    decimal_precision: int,
    max_steps: int = 10000,
    extra_precision: int = 200,
) -> list[mp.mpc]:
    """Return every finite root using only arbitrary-precision arithmetic."""

    mp.mp.dps = decimal_precision
    retained = list(coefficients)
    while len(retained) > 1 and retained[-1] == 0:
        retained.pop()
    if len(retained) <= 1:
        return []

    # Exact zero constant terms encode roots at s=0.  Deflate them before
    # calling Durand--Kerner: mp.polyroots converges poorly for a high-
    # multiplicity zero (notably the exact AG=0 row), while the deflation is
    # algebraically exact and retains the full root multiplicity.
    zero_multiplicity = 0
    while len(retained) > 1 and retained[0] == 0:
        retained.pop(0)
        zero_multiplicity += 1

    # Balance in s = alpha*z using arbitrary-precision coefficient ratios.
    degree = len(retained) - 1
    leading = abs(retained[-1])
    estimates = [
        (abs(value) / leading) ** (mp.mpf(1) / (degree - power))
        for power, value in enumerate(retained[:-1])
        if value != 0
    ]
    alpha = (
        mp.exp(sum(mp.log(value) for value in estimates) / len(estimates))
        if estimates else mp.mpf(1)
    )
    scaled = [value * alpha**power for power, value in enumerate(retained)]
    normalization = max(abs(value) for value in scaled)
    roots_z = mp.polyroots(
        [value / normalization for value in reversed(scaled)],
        maxsteps=max_steps,
        error=False,
        extraprec=extra_precision,
    )
    return [mp.mpc(0)] * zero_multiplicity + [alpha * root for root in roots_z]


def _mp_poly_value(s: mp.mpc, coefficients: list[mp.mpc]) -> mp.mpc:
    value = mp.mpc(0)
    for coefficient in reversed(coefficients):
        value = value * s + coefficient
    return value


def _mp_poly_derivative_value(s: mp.mpc, coefficients: list[mp.mpc]) -> mp.mpc:
    value = mp.mpc(0)
    for power in range(len(coefficients) - 1, 0, -1):
        value = value * s + power * coefficients[power]
    return value


def _mp_dstar_and_scale(
    s: mp.mpc,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
) -> tuple[mp.mpc, mp.mpf]:
    arho = mp.mpf(str(float(Arho)))
    amu = mp.mpf(str(float(Amu)))
    ag = mp.mpf(str(float(AG)))
    lam = mp.mpf(str(float(Lambda)))

    eta1 = (1 + amu) + (lam / s) * (1 + ag)
    eta2 = (1 - amu) + (lam / s) * (1 - ag)
    q1 = mp.sqrt(1 + ((1 + arho) * s) / eta1)
    q2 = mp.sqrt(1 + ((1 - arho) * s) / eta2)
    denominator1 = eta1 + eta2 * q2
    denominator2 = eta2 + eta1 * q1
    term1 = s / denominator1
    term2 = s / denominator2
    value = term1 + term2 + 2
    scale = abs(term1) + abs(term2) + 2
    return value, scale


def _threshold(
    residual_scale: mp.mpf,
    *,
    relative_tolerance: float,
    absolute_tolerance: float,
) -> mp.mpf:
    return mp.mpf(str(float(absolute_tolerance))) + mp.mpf(
        str(float(relative_tolerance))
    ) * max(mp.mpf(1), residual_scale)


def _nearest_spacing(candidate: complex, candidates: Iterable[complex]) -> float | None:
    distances = [
        abs(complex(candidate) - complex(other))
        for other in candidates
        if complex(candidate) != complex(other)
    ]
    if not distances:
        return None
    return float(min(distances))


def _cluster_size(
    candidate: complex,
    candidates: Iterable[complex],
    *,
    relative_threshold: float,
) -> int:
    candidate = complex(candidate)
    count = 0
    for other in candidates:
        other = complex(other)
        threshold = relative_threshold * max(1.0, abs(candidate), abs(other))
        if abs(candidate - other) <= threshold:
            count += 1
    return count


def _newton_refine(
    initial: complex,
    coefficients: list[mp.mpc],
    *,
    max_iterations: int,
) -> mp.mpc:
    root = mp.mpc(str(float(complex(initial).real)), str(float(complex(initial).imag)))
    for _ in range(max_iterations):
        value = _mp_poly_value(root, coefficients)
        derivative = _mp_poly_derivative_value(root, coefficients)
        if abs(derivative) == 0:
            raise ZeroDivisionError("zero derivative during high-precision polishing")
        step = value / derivative
        root -= step
        if abs(step) <= mp.eps * max(mp.mpf(1), abs(root)):
            break
    return root


def _full_solve_nearest(
    initial: complex,
    coefficients: list[mp.mpc],
    root_scale: float,
    *,
    max_steps: int,
) -> mp.mpc:
    alpha = mp.mpf(str(float(root_scale)))
    powers = [alpha**power for power in range(len(coefficients))]
    scaled_ascending = [
        coefficient * powers[power]
        for power, coefficient in enumerate(coefficients)
    ]
    roots_z = mp.polyroots(
        list(reversed(scaled_ascending)),
        maxsteps=max_steps,
        error=False,
    )
    target_z = mp.mpc(str(float(complex(initial).real)), str(float(complex(initial).imag))) / alpha
    nearest_z = min(roots_z, key=lambda root: abs(root - target_z))
    return alpha * nearest_z


def refine_ambiguous_candidate(
    verification_result: CandidateVerificationResult,
    all_candidates: Iterable[complex],
    effective_polynomial: EffectivePolynomial,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    residual_relative_tolerance: float,
    residual_absolute_tolerance: float,
    ambiguity_factor: float,
    options: HighPrecisionOptions = HighPrecisionOptions(),
    domain_relative_tolerance: float = 1e-10,
    domain_absolute_tolerance: float = 1e-14,
) -> tuple[CandidateVerificationResult, HighPrecisionRefinementResult]:
    """Refine one ambiguous candidate and return an updated verification row."""

    original = complex(verification_result.candidate)
    nearest_spacing = _nearest_spacing(original, all_candidates)
    cluster_size = _cluster_size(
        original,
        all_candidates,
        relative_threshold=options.cluster_spacing_threshold,
    )

    last_reason = "not attempted"
    precision = options.initial_decimal_precision
    while precision <= options.maximum_decimal_precision:
        try:
            mp.mp.dps = precision
            coefficients = mp_p14_coefficients(
                Arho,
                Amu,
                AG,
                Lambda,
                decimal_precision=precision,
            )
            full_solve_used = False
            try:
                refined = _newton_refine(
                    original,
                    coefficients,
                    max_iterations=options.maximum_refinement_attempts,
                )
            except Exception:
                refined = _full_solve_nearest(
                    original,
                    coefficients,
                    effective_polynomial.root_scale,
                    max_steps=options.maximum_refinement_attempts,
                )
                full_solve_used = True
            p14_residual = abs(_mp_poly_value(refined, coefficients))
            dstar_value, residual_scale = _mp_dstar_and_scale(
                refined,
                Arho,
                Amu,
                AG,
                Lambda,
            )
            absolute_residual = abs(dstar_value)
            relative_residual = absolute_residual / max(mp.mpf(1), residual_scale)
            threshold = _threshold(
                residual_scale,
                relative_tolerance=residual_relative_tolerance,
                absolute_tolerance=residual_absolute_tolerance,
            )

            refined_complex = complex(float(mp.re(refined)), float(mp.im(refined)))
            domain_result = check_candidate_domain(
                refined_complex,
                Arho,
                Amu,
                AG,
                Lambda,
                relative_tolerance=domain_relative_tolerance,
                absolute_tolerance=domain_absolute_tolerance,
            )

            if not domain_result.domain_valid:
                status: RefinementStatus = "still_ambiguous"
                reason = "refined_candidate_failed_domain_check"
                candidate_status = "domain_invalid"
                mathematically_genuine = False
                rejection_reason = "domain failure after high-precision refinement"
            elif absolute_residual <= threshold:
                status = "refined_genuine"
                reason = "high_precision_residual_satisfied"
                candidate_status = "mathematically_genuine"
                mathematically_genuine = True
                rejection_reason = None
            elif absolute_residual <= mp.mpf(str(float(ambiguity_factor))) * threshold:
                status = "still_ambiguous"
                reason = "high_precision_residual_near_threshold"
                candidate_status = "numerically_ambiguous"
                mathematically_genuine = False
                rejection_reason = "residual near verification threshold after refinement"
            else:
                status = "refined_spurious"
                reason = "high_precision_residual_too_large"
                candidate_status = "residual_too_large"
                mathematically_genuine = False
                rejection_reason = "residual too large after high-precision refinement"

            hp_result = HighPrecisionRefinementResult(
                original_candidate=original,
                refined_candidate=refined_complex,
                working_precision=precision,
                high_precision_p14_residual=float(p14_residual),
                high_precision_absolute_residual=float(absolute_residual),
                high_precision_relative_residual=float(relative_residual),
                verification_threshold=float(threshold),
                status=status,
                reason=reason,
                cluster_size=cluster_size,
                nearest_root_spacing=nearest_spacing,
                full_polynomial_solve_used=full_solve_used,
            )
            updated = CandidateVerificationResult(
                candidate=refined_complex,
                domain_result=domain_result,
                domain_valid=domain_result.domain_valid,
                absolute_residual=float(absolute_residual),
                relative_residual=float(relative_residual),
                residual_scale=float(residual_scale),
                verification_threshold=float(threshold),
                mathematically_genuine=mathematically_genuine,
                status=candidate_status,  # type: ignore[arg-type]
                rejection_reason=None if mathematically_genuine else rejection_reason,
                dispersion_value=complex(float(mp.re(dstar_value)), float(mp.im(dstar_value))),
                high_precision_refinement=hp_result,
            )
            return updated, hp_result
        except Exception as exc:  # noqa: BLE001 - preserve refinement diagnostics
            last_reason = f"{type(exc).__name__}: {exc}"
            precision += options.precision_step

    hp_result = HighPrecisionRefinementResult(
        original_candidate=original,
        refined_candidate=None,
        working_precision=min(precision, options.maximum_decimal_precision),
        high_precision_p14_residual=None,
        high_precision_absolute_residual=None,
        high_precision_relative_residual=None,
        verification_threshold=None,
        status="refinement_failed",
        reason=last_reason,
        cluster_size=cluster_size,
        nearest_root_spacing=nearest_spacing,
        full_polynomial_solve_used=False,
    )
    return (
        CandidateVerificationResult(
            candidate=verification_result.candidate,
            domain_result=verification_result.domain_result,
            domain_valid=verification_result.domain_valid,
            absolute_residual=verification_result.absolute_residual,
            relative_residual=verification_result.relative_residual,
            residual_scale=verification_result.residual_scale,
            verification_threshold=verification_result.verification_threshold,
            mathematically_genuine=False,
            status="numerically_ambiguous",
            rejection_reason="high-precision refinement failed",
            dispersion_value=verification_result.dispersion_value,
            high_precision_refinement=hp_result,
        ),
        hp_result,
    )


def refine_ambiguous_verification_run(
    verification_run: VerificationRun,
    effective_polynomial: EffectivePolynomial,
    Arho: float,
    Amu: float,
    AG: float,
    Lambda: float,
    *,
    residual_relative_tolerance: float,
    residual_absolute_tolerance: float,
    ambiguity_factor: float,
    options: HighPrecisionOptions = HighPrecisionOptions(),
    domain_relative_tolerance: float = 1e-10,
    domain_absolute_tolerance: float = 1e-14,
) -> VerificationRun:
    """Resolve ambiguous records and genuine members of their root clusters.

    Root estimates for a tightly clustered polynomial can move enough between
    LAPACK/NumPy versions for one member to fall just above or below the
    double-precision verification threshold.  Once a cluster contains an
    ambiguous member, refine its provisionally genuine members as well so the
    cluster is treated consistently.  Clearly spurious members are not sent
    through the expensive fallback.
    """

    candidates = tuple(result.candidate for result in verification_run.results)
    ambiguous_candidates = tuple(
        result.candidate
        for result in verification_run.results
        if result.status == "numerically_ambiguous"
    )

    def belongs_to_ambiguous_cluster(result: CandidateVerificationResult) -> bool:
        if result.status == "numerically_ambiguous":
            return True
        if result.status != "mathematically_genuine":
            return False
        candidate = complex(result.candidate)
        return any(
            abs(candidate - complex(ambiguous))
            <= options.cluster_spacing_threshold
            * max(1.0, abs(candidate), abs(complex(ambiguous)))
            for ambiguous in ambiguous_candidates
        )

    refined_results: list[CandidateVerificationResult] = []
    for result in verification_run.results:
        if not belongs_to_ambiguous_cluster(result):
            refined_results.append(result)
            continue
        refined, _diagnostic = refine_ambiguous_candidate(
            result,
            candidates,
            effective_polynomial,
            Arho,
            Amu,
            AG,
            Lambda,
            residual_relative_tolerance=residual_relative_tolerance,
            residual_absolute_tolerance=residual_absolute_tolerance,
            ambiguity_factor=ambiguity_factor,
            options=options,
            domain_relative_tolerance=domain_relative_tolerance,
            domain_absolute_tolerance=domain_absolute_tolerance,
        )
        refined_results.append(refined)

    result_tuple = tuple(refined_results)
    return VerificationRun(
        results=result_tuple,
        summary=summarize_verification(result_tuple),
    )
