"""Batch parameter-sweep infrastructure.

This module evaluates collections of nondimensional parameter tuples by calling
``src.workflow.analyze_parameter_point`` as the atomic operation.  It does not
reimplement polynomial generation, domain checks, mathematical verification,
or physical admissibility.

Checkpoint format
-----------------

Checkpoints are newline-delimited JSON (JSON Lines).  Each line is either a
``metadata`` record or a completed ``point`` record.  Point records contain
JSON-safe summaries plus serialized complex roots and diagnostics.  The in-memory
``SweepResult`` preserves the full ``ParameterPointResult`` object for every
point analyzed in the current run.
"""

from __future__ import annotations

from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
import itertools
import json
from pathlib import Path
import time
from typing import Any, Callable, Iterable, Literal, Mapping

import numpy as np

from .workflow import (
    ParameterClassification,
    ParameterPointResult,
    analyze_parameter_point,
)


WORKFLOW_VERSION = "p14-root-workflow-sprint6-high-precision-refinement"
CLASSIFICATIONS: tuple[ParameterClassification, ...] = (
    "outside_theorem_scope",
    "identically_zero_polynomial",
    "constant_nonzero_polynomial",
    "only_domain_invalid_candidates",
    "only_spurious_candidates",
    "genuine_roots_nonadmissible",
    "marginal_spatial_decay",
    "physically_admissible_roots",
    "numerically_unresolved",
)


@dataclass(frozen=True)
class SweepPointInput:
    """One user-facing point in a sweep definition."""

    index: int
    coordinates: dict[str, float]
    parameters: dict[str, float]


@dataclass(frozen=True)
class SweepDefinition:
    """Machine-readable description of the requested sweep."""

    kind: Literal["explicit", "cartesian_grid"]
    points: tuple[SweepPointInput, ...]
    fixed_values: dict[str, float]
    swept_values: dict[str, tuple[float, ...]]
    coordinate_system: Literal["Lambda", "Ek", "mixed"]


@dataclass(frozen=True)
class SweepPointResult:
    """Result for one sweep point."""

    index: int
    coordinates: dict[str, float]
    parameters: dict[str, float]
    result: ParameterPointResult | None
    classification: ParameterClassification
    runtime_seconds: float
    error: str | None = None
    traceback_summary: str | None = None
    skipped_existing: bool = False
    checkpoint_record: dict[str, Any] | None = None


@dataclass(frozen=True)
class SweepResult:
    """Complete batch sweep result."""

    definition: SweepDefinition
    points: tuple[SweepPointResult, ...]
    aggregate_classification_counts: dict[str, int]
    unresolved_points: tuple[int, ...]
    total_runtime_seconds: float
    per_point_runtime_seconds: tuple[float, ...]
    checkpoint_path: str | None
    workflow_version: str = WORKFLOW_VERSION


def lambda_from_ek(Ek: float) -> float:
    """Convert bounded sweep coordinate ``Ek`` to canonical ``Lambda``.

    The accepted finite coordinate range is ``-1 < Ek <= 1``.  The endpoint
    ``Ek = 1`` maps to ``Lambda = 0``, which is outside the current no-root-loss
    theorem scope but is still represented explicitly.
    """

    ek = float(Ek)
    if not np.isfinite(ek):
        raise ValueError("Ek must be finite")
    if not (-1.0 < ek <= 1.0):
        raise ValueError("Ek must lie in (-1, 1]")
    return float(((1.0 - ek) / (1.0 + ek)) ** 2)


def _canonicalize_point(index: int, point: Mapping[str, float]) -> SweepPointInput:
    coordinates = {key: float(value) for key, value in point.items()}
    has_lambda = "Lambda" in coordinates
    has_ek = "Ek" in coordinates
    if has_lambda and has_ek:
        converted = lambda_from_ek(coordinates["Ek"])
        if not np.isclose(converted, coordinates["Lambda"], rtol=1e-12, atol=1e-14):
            raise ValueError("Provided Ek and Lambda are inconsistent")
        lam = float(coordinates["Lambda"])
    elif has_ek:
        lam = lambda_from_ek(coordinates["Ek"])
        coordinates["Lambda"] = lam
    elif has_lambda:
        lam = float(coordinates["Lambda"])
    else:
        raise ValueError("Each point must provide either Lambda or Ek")

    required = ("Arho", "Amu", "AG")
    missing = [name for name in required if name not in coordinates]
    if missing:
        raise ValueError("Missing required parameter(s): " + ", ".join(missing))

    parameters = {
        "Arho": float(coordinates["Arho"]),
        "Amu": float(coordinates["Amu"]),
        "AG": float(coordinates["AG"]),
        "Lambda": float(lam),
    }
    return SweepPointInput(
        index=index,
        coordinates=coordinates,
        parameters=parameters,
    )


def explicit_sweep_definition(
    points: Iterable[Mapping[str, float]],
) -> SweepDefinition:
    """Build a sweep definition from an explicit list of point dictionaries."""

    sweep_points = tuple(
        _canonicalize_point(index, point)
        for index, point in enumerate(points)
    )
    coordinate_system = _coordinate_system_for_points(sweep_points)
    return SweepDefinition(
        kind="explicit",
        points=sweep_points,
        fixed_values={},
        swept_values={},
        coordinate_system=coordinate_system,
    )


def cartesian_grid_definition(
    *,
    fixed: Mapping[str, float],
    sweep: Mapping[str, Iterable[float]],
) -> SweepDefinition:
    """Build a Cartesian-grid sweep definition."""

    sweep_values = {
        key: tuple(float(value) for value in values)
        for key, values in sweep.items()
    }
    for key, values in sweep_values.items():
        if not values:
            raise ValueError(f"sweep values for {key} must be nonempty")

    keys = tuple(sweep_values)
    raw_points = []
    for combination in itertools.product(*(sweep_values[key] for key in keys)):
        point = {key: float(value) for key, value in fixed.items()}
        point.update(dict(zip(keys, combination, strict=True)))
        raw_points.append(point)

    sweep_points = tuple(
        _canonicalize_point(index, point)
        for index, point in enumerate(raw_points)
    )
    return SweepDefinition(
        kind="cartesian_grid",
        points=sweep_points,
        fixed_values={key: float(value) for key, value in fixed.items()},
        swept_values=sweep_values,
        coordinate_system=_coordinate_system_for_points(sweep_points),
    )


def _coordinate_system_for_points(
    points: tuple[SweepPointInput, ...],
) -> Literal["Lambda", "Ek", "mixed"]:
    has_ek = any("Ek" in point.coordinates for point in points)
    has_lambda_without_ek = any(
        "Lambda" in point.coordinates and "Ek" not in point.coordinates
        for point in points
    )
    if has_ek and has_lambda_without_ek:
        return "mixed"
    if has_ek:
        return "Ek"
    return "Lambda"


def _complex_to_json(value: complex | None) -> dict[str, float] | None:
    if value is None:
        return None
    z = complex(value)
    return {"real": float(z.real), "imag": float(z.imag)}


def _complex_list(values: Iterable[complex]) -> list[dict[str, float]]:
    return [_complex_to_json(value) for value in values]  # type: ignore[list-item]


def _counts_to_json(result: ParameterPointResult) -> dict[str, Any]:
    counts = result.counts
    return {
        "effective_polynomial_degree": counts.effective_polynomial_degree,
        "total_polynomial_candidates": counts.total_polynomial_candidates,
        "domain_valid_candidates": counts.domain_valid_candidates,
        "domain_invalid_candidates": counts.domain_invalid_candidates,
        "mathematically_genuine_roots": counts.mathematically_genuine_roots,
        "ambiguous_mathematical_candidates": counts.ambiguous_mathematical_candidates,
        "physically_admissible_roots": counts.physically_admissible_roots,
        "marginal_spatial_decay_roots": counts.marginal_spatial_decay_roots,
        "physically_nonadmissible_roots": counts.physically_nonadmissible_roots,
        "evaluation_failures": counts.evaluation_failures,
    }


def parameter_point_to_json(result: ParameterPointResult) -> dict[str, Any]:
    """Serialize a single-point result to JSON-compatible data."""

    effective = result.effective_polynomial
    domain_rows = []
    for row in result.domain_results:
        intermediates = row.intermediates
        domain_rows.append(
            {
                "candidate": _complex_to_json(row.s),
                "domain_valid": row.domain_valid,
                "failed_conditions": list(row.failed_conditions),
                "notes": list(row.notes),
                "intermediates": None
                if intermediates is None
                else {
                    "E_plus": _complex_to_json(intermediates.E_plus),
                    "E_minus": _complex_to_json(intermediates.E_minus),
                    "R_plus": _complex_to_json(intermediates.R_plus),
                    "R_minus": _complex_to_json(intermediates.R_minus),
                    "q_plus": _complex_to_json(intermediates.q_plus),
                    "q_minus": _complex_to_json(intermediates.q_minus),
                    "B_plus": _complex_to_json(intermediates.B_plus),
                    "B_minus": _complex_to_json(intermediates.B_minus),
                },
                "checks": [
                    {
                        "name": check.name,
                        "value": _complex_to_json(check.value),
                        "magnitude": check.magnitude,
                        "threshold": check.threshold,
                        "passed": check.passed,
                    }
                    for check in row.checks
                ],
            }
        )

    verification_rows = []
    if result.verification_run is not None:
        for row in result.verification_run.results:
            refinement = row.high_precision_refinement
            verification_rows.append(
                {
                    "candidate": _complex_to_json(row.candidate),
                    "domain_valid": row.domain_valid,
                    "absolute_residual": row.absolute_residual,
                    "relative_residual": row.relative_residual,
                    "residual_scale": row.residual_scale,
                    "verification_threshold": row.verification_threshold,
                    "mathematically_genuine": row.mathematically_genuine,
                    "status": row.status,
                    "rejection_reason": row.rejection_reason,
                    "dispersion_value": _complex_to_json(row.dispersion_value),
                    "high_precision_refinement": None
                    if refinement is None
                    else {
                        "original_candidate": _complex_to_json(
                            refinement.original_candidate
                        ),
                        "refined_candidate": _complex_to_json(
                            refinement.refined_candidate
                        ),
                        "working_precision": refinement.working_precision,
                        "high_precision_p14_residual": (
                            refinement.high_precision_p14_residual
                        ),
                        "high_precision_absolute_residual": (
                            refinement.high_precision_absolute_residual
                        ),
                        "high_precision_relative_residual": (
                            refinement.high_precision_relative_residual
                        ),
                        "verification_threshold": refinement.verification_threshold,
                        "status": refinement.status,
                        "reason": refinement.reason,
                        "cluster_size": refinement.cluster_size,
                        "nearest_root_spacing": refinement.nearest_root_spacing,
                        "full_polynomial_solve_used": (
                            refinement.full_polynomial_solve_used
                        ),
                    },
                }
            )

    physical_rows = []
    if result.physical_run is not None:
        for row in result.physical_run.results:
            physical_rows.append(
                {
                    "candidate": _complex_to_json(row.candidate),
                    "mathematically_genuine": row.mathematically_genuine,
                    "q_plus": _complex_to_json(row.q_plus),
                    "q_minus": _complex_to_json(row.q_minus),
                    "re_q_plus": row.re_q_plus,
                    "re_q_minus": row.re_q_minus,
                    "tolerance_plus": row.tolerance_plus,
                    "tolerance_minus": row.tolerance_minus,
                    "physically_admissible": row.physically_admissible,
                    "status": row.status,
                    "reason": row.reason,
                }
            )

    return {
        "parameters": {
            "Arho": result.parameters.Arho,
            "Amu": result.parameters.Amu,
            "AG": result.parameters.AG,
            "Lambda": result.parameters.Lambda,
        },
        "theorem_scope": {
            "applies": result.theorem_scope.applies,
            "failed_conditions": list(result.theorem_scope.failed_conditions),
            "analysis_skipped": result.theorem_scope.analysis_skipped,
            "diagnostic_only": result.theorem_scope.diagnostic_only,
        },
        "effective_polynomial": None
        if effective is None
        else {
            "status": effective.status,
            "effective_degree": effective.effective_degree,
            "coefficient_scale": effective.coefficient_scale,
            "trim_tolerance": effective.trim_tolerance,
            "normalization_factor": effective.normalization_factor,
            "root_scale": effective.root_scale,
            "raw_coefficients": _complex_list(effective.raw_coefficients),
        },
        "reduced_degree": result.reduced_degree,
        "classification": result.classification,
        "warnings": list(result.warnings),
        "counts": _counts_to_json(result),
        "polynomial_candidates": _complex_list(result.polynomial_candidates),
        "domain_results": domain_rows,
        "verification_results": verification_rows,
        "physical_results": physical_rows,
        "physical_root_metadata": [
            {
                "root": _complex_to_json(row.root),
                "temporal_status": row.temporal_status,
                "real_part": row.real_part,
                "imaginary_part": row.imaginary_part,
            }
            for row in result.physical_root_metadata
        ],
    }


def sweep_point_to_json(point: SweepPointResult) -> dict[str, Any]:
    """Serialize one sweep point to a JSON-compatible record."""

    return {
        "record_type": "point",
        "workflow_version": WORKFLOW_VERSION,
        "index": point.index,
        "coordinates": point.coordinates,
        "parameters": point.parameters,
        "classification": point.classification,
        "runtime_seconds": point.runtime_seconds,
        "error": point.error,
        "traceback_summary": point.traceback_summary,
        "skipped_existing": point.skipped_existing,
        "checkpoint_record": point.checkpoint_record,
        "result": None if point.result is None else parameter_point_to_json(point.result),
    }


def sweep_definition_to_json(definition: SweepDefinition) -> dict[str, Any]:
    """Serialize a sweep definition to JSON-compatible metadata."""

    return {
        "kind": definition.kind,
        "fixed_values": definition.fixed_values,
        "swept_values": {
            key: list(values) for key, values in definition.swept_values.items()
        },
        "coordinate_system": definition.coordinate_system,
        "point_count": len(definition.points),
        "points": [
            {
                "index": point.index,
                "coordinates": point.coordinates,
                "parameters": point.parameters,
            }
            for point in definition.points
        ],
    }


def _checkpoint_metadata(definition: SweepDefinition, analyze_kwargs: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "record_type": "metadata",
        "workflow_version": WORKFLOW_VERSION,
        "created_unix_time": time.time(),
        "definition": sweep_definition_to_json(definition),
        "analyze_kwargs": dict(analyze_kwargs),
    }


def _append_jsonl(path: Path, record: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def load_checkpoint_records(path: str | Path) -> list[dict[str, Any]]:
    """Load JSONL checkpoint records."""

    checkpoint = Path(path)
    if not checkpoint.exists():
        return []
    records = []
    with checkpoint.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def completed_indices_from_checkpoint(path: str | Path) -> set[int]:
    """Return sweep-point indices already present in a checkpoint."""

    return {
        int(record["index"])
        for record in load_checkpoint_records(path)
        if record.get("record_type") == "point"
    }


def completed_point_records_from_checkpoint(path: str | Path) -> dict[int, dict[str, Any]]:
    """Return latest checkpoint point record for each completed index."""

    records: dict[int, dict[str, Any]] = {}
    for record in load_checkpoint_records(path):
        if record.get("record_type") == "point":
            records[int(record["index"])] = record
    return records


def _analyze_one(
    point: SweepPointInput,
    analyze_kwargs: Mapping[str, Any],
) -> SweepPointResult:
    start = time.perf_counter()
    try:
        result = analyze_parameter_point(
            point.parameters["Arho"],
            point.parameters["Amu"],
            point.parameters["AG"],
            point.parameters["Lambda"],
            **dict(analyze_kwargs),
        )
        return SweepPointResult(
            index=point.index,
            coordinates=point.coordinates,
            parameters=point.parameters,
            result=result,
            classification=result.classification,
            runtime_seconds=time.perf_counter() - start,
        )
    except Exception as exc:  # noqa: BLE001 - point failures must not kill sweep
        return SweepPointResult(
            index=point.index,
            coordinates=point.coordinates,
            parameters=point.parameters,
            result=None,
            classification="numerically_unresolved",
            runtime_seconds=time.perf_counter() - start,
            error=type(exc).__name__,
            traceback_summary=str(exc),
        )


def _worker_analyze_one(
    point: SweepPointInput,
    analyze_kwargs: Mapping[str, Any],
) -> SweepPointResult:
    return _analyze_one(point, analyze_kwargs)


def aggregate_classification_counts(
    points: Iterable[SweepPointResult],
) -> dict[str, int]:
    """Aggregate mutually exclusive final classifications."""

    counter = Counter(point.classification for point in points)
    return {classification: int(counter.get(classification, 0)) for classification in CLASSIFICATIONS}


def run_parameter_sweep(
    definition: SweepDefinition,
    *,
    max_workers: int = 1,
    checkpoint_path: str | Path | None = None,
    resume: bool = False,
    progress_callback: Callable[[SweepPointResult, dict[str, int]], None] | None = None,
    **analyze_kwargs: Any,
) -> SweepResult:
    """Evaluate every point in a sweep definition.

    Serial execution is the default.  If ``max_workers > 1``, a process pool is
    used.  Final in-memory results are sorted by input index regardless of
    execution order.
    """

    if max_workers < 1:
        raise ValueError("max_workers must be at least one")

    start = time.perf_counter()
    checkpoint = None if checkpoint_path is None else Path(checkpoint_path)
    completed: set[int] = set()
    skipped_results: list[SweepPointResult] = []
    if checkpoint is not None:
        if resume:
            completed_records = completed_point_records_from_checkpoint(checkpoint)
            completed = set(completed_records)
            for index, record in completed_records.items():
                skipped_results.append(
                    SweepPointResult(
                        index=index,
                        coordinates=dict(record.get("coordinates", {})),
                        parameters=dict(record.get("parameters", {})),
                        result=None,
                        classification=record.get("classification", "numerically_unresolved"),
                        runtime_seconds=0.0,
                        error=record.get("error"),
                        traceback_summary=record.get("traceback_summary"),
                        skipped_existing=True,
                        checkpoint_record=record,
                    )
                )
        else:
            checkpoint.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.write_text("", encoding="utf-8")
        _append_jsonl(checkpoint, _checkpoint_metadata(definition, analyze_kwargs))

    pending = [point for point in definition.points if point.index not in completed]
    results: list[SweepPointResult] = list(skipped_results)

    if max_workers == 1:
        for point in pending:
            result = _analyze_one(point, analyze_kwargs)
            results.append(result)
            if checkpoint is not None:
                _append_jsonl(checkpoint, sweep_point_to_json(result))
            if progress_callback is not None:
                progress_callback(result, aggregate_classification_counts(results))
    else:
        try:
            executor_context = ProcessPoolExecutor(max_workers=max_workers)
        except (OSError, PermissionError):
            executor_context = None

        if executor_context is None:
            for point in pending:
                result = _analyze_one(point, analyze_kwargs)
                results.append(result)
                if checkpoint is not None:
                    _append_jsonl(checkpoint, sweep_point_to_json(result))
                if progress_callback is not None:
                    progress_callback(result, aggregate_classification_counts(results))
        else:
            with executor_context as executor:
                futures = {
                    executor.submit(_worker_analyze_one, point, analyze_kwargs): point.index
                    for point in pending
                }
                for future in as_completed(futures):
                    result = future.result()
                    results.append(result)
                    if checkpoint is not None:
                        _append_jsonl(checkpoint, sweep_point_to_json(result))
                    if progress_callback is not None:
                        progress_callback(result, aggregate_classification_counts(results))

    results.sort(key=lambda row: row.index)
    counts = aggregate_classification_counts(results)
    unresolved = tuple(
        point.index
        for point in results
        if point.classification == "numerically_unresolved"
    )
    runtimes = tuple(point.runtime_seconds for point in results)
    return SweepResult(
        definition=definition,
        points=tuple(results),
        aggregate_classification_counts=counts,
        unresolved_points=unresolved,
        total_runtime_seconds=time.perf_counter() - start,
        per_point_runtime_seconds=runtimes,
        checkpoint_path=None if checkpoint is None else str(checkpoint),
    )
