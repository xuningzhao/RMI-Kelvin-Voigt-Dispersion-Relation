"""Postprocess Sprint 6 parameter-sweep JSONL outputs.

This module is intentionally read-only with respect to the mathematical
workflow: it never calls ``analyze_parameter_point`` and never recomputes
polynomial roots, dispersion residuals, domain checks, or physical
admissibility.  It only extracts compact, visualization-ready tables from the
full JSONL diagnostics produced by ``src.parameter_sweep``.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import csv
import itertools
import json
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

from .parameter_sweep import CLASSIFICATIONS


CLASSIFICATION_CODES: dict[str, int] = {
    classification: index for index, classification in enumerate(CLASSIFICATIONS)
}
"""Stable categorical encoding for final parameter-point classifications.

The codes are identifiers only.  They do not imply severity, stability,
ordering, or physical preference.
"""


POINT_SUMMARY_FIELDS = (
    "point_id",
    "Arho",
    "Amu",
    "AG",
    "Lambda",
    "Ek",
    "theorem_scope_applies",
    "theorem_scope_failed_conditions",
    "final_classification",
    "classification_code",
    "effective_polynomial_degree",
    "reduced_degree",
    "polynomial_candidate_count",
    "domain_valid_candidate_count",
    "mathematically_genuine_root_count",
    "physically_admissible_root_count",
    "marginal_root_count",
    "physically_nonadmissible_root_count",
    "numerically_unresolved",
    "runtime_seconds",
    "warning_count",
)


ROOT_SUMMARY_FIELDS = (
    "point_id",
    "root_index",
    "s_real",
    "s_imag",
    "temporal_status",
    "q_plus_real",
    "q_plus_imag",
    "q_minus_real",
    "q_minus_imag",
    "absolute_residual",
    "relative_residual",
    "high_precision_refinement_status",
    "admissibility_status",
)


class SweepSummaryError(ValueError):
    """Raised when a Sprint 6 JSONL file cannot be summarized safely."""


@dataclass(frozen=True)
class SweepSummary:
    """Compact representation extracted from a Sprint 6 JSONL file."""

    source_jsonl_path: str
    workflow_version: str | None
    sweep_definition: dict[str, Any]
    analysis_options: dict[str, Any]
    creation_time: float | None
    parameter_ordering: tuple[str, ...]
    coordinate_system: str
    point_rows: tuple[dict[str, Any], ...]
    admissible_root_rows: tuple[dict[str, Any], ...]
    aggregate_classification_counts: dict[str, int]
    validation_report: dict[str, Any]
    metadata_records: tuple[dict[str, Any], ...]


def _load_jsonl(path: str | Path) -> list[dict[str, Any]]:
    source = Path(path)
    if not source.exists():
        raise SweepSummaryError(f"Input JSONL file does not exist: {source}")
    records: list[dict[str, Any]] = []
    with source.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped:
                continue
            try:
                record = json.loads(stripped)
            except json.JSONDecodeError as exc:
                raise SweepSummaryError(
                    f"Malformed JSON on line {line_number}: {exc}"
                ) from exc
            if not isinstance(record, dict):
                raise SweepSummaryError(f"Line {line_number} is not a JSON object")
            records.append(record)
    if not records:
        raise SweepSummaryError("Input JSONL file is empty")
    return records


def _complex_parts(value: Mapping[str, Any] | None) -> tuple[float | None, float | None]:
    if value is None:
        return None, None
    return float(value["real"]), float(value["imag"])


def _parameter_value(
    coordinates: Mapping[str, Any],
    parameters: Mapping[str, Any],
    name: str,
) -> float | None:
    if name in coordinates:
        return float(coordinates[name])
    if name in parameters:
        return float(parameters[name])
    return None


def _semicolon_join(values: Iterable[Any]) -> str:
    return ";".join(str(value) for value in values)


def _extract_point_row(record: Mapping[str, Any]) -> dict[str, Any]:
    result = record.get("result") or {}
    counts = result.get("counts") or {}
    theorem_scope = result.get("theorem_scope") or {}
    coordinates = record.get("coordinates") or result.get("coordinates") or {}
    parameters = record.get("parameters") or result.get("parameters") or {}
    classification = str(record.get("classification") or result.get("classification"))

    if classification not in CLASSIFICATION_CODES:
        raise SweepSummaryError(f"Unknown classification: {classification}")

    warnings = result.get("warnings") or []
    return {
        "point_id": int(record["index"]),
        "Arho": _parameter_value(coordinates, parameters, "Arho"),
        "Amu": _parameter_value(coordinates, parameters, "Amu"),
        "AG": _parameter_value(coordinates, parameters, "AG"),
        "Lambda": _parameter_value(coordinates, parameters, "Lambda"),
        "Ek": _parameter_value(coordinates, parameters, "Ek"),
        "theorem_scope_applies": bool(theorem_scope.get("applies", False)),
        "theorem_scope_failed_conditions": _semicolon_join(
            theorem_scope.get("failed_conditions") or []
        ),
        "final_classification": classification,
        "classification_code": CLASSIFICATION_CODES[classification],
        "effective_polynomial_degree": counts.get("effective_polynomial_degree"),
        "reduced_degree": bool(result.get("reduced_degree", False)),
        "polynomial_candidate_count": int(counts.get("total_polynomial_candidates", 0)),
        "domain_valid_candidate_count": int(counts.get("domain_valid_candidates", 0)),
        "mathematically_genuine_root_count": int(
            counts.get("mathematically_genuine_roots", 0)
        ),
        "physically_admissible_root_count": int(
            counts.get("physically_admissible_roots", 0)
        ),
        "marginal_root_count": int(counts.get("marginal_spatial_decay_roots", 0)),
        "physically_nonadmissible_root_count": int(
            counts.get("physically_nonadmissible_roots", 0)
        ),
        "numerically_unresolved": classification == "numerically_unresolved",
        "runtime_seconds": float(record.get("runtime_seconds", 0.0)),
        "warning_count": len(warnings),
    }


def _extract_root_rows(record: Mapping[str, Any]) -> list[dict[str, Any]]:
    result = record.get("result") or {}
    physical_rows = result.get("physical_results") or []
    verification_rows = result.get("verification_results") or []
    rows: list[dict[str, Any]] = []
    for local_index, physical in enumerate(physical_rows):
        if not physical.get("physically_admissible", False):
            continue
        s_real, s_imag = _complex_parts(physical.get("candidate"))
        q_plus_real, q_plus_imag = _complex_parts(physical.get("q_plus"))
        q_minus_real, q_minus_imag = _complex_parts(physical.get("q_minus"))
        verification = (
            verification_rows[local_index]
            if local_index < len(verification_rows)
            else {}
        )
        rows.append(
            {
                "point_id": int(record["index"]),
                "root_index": len(rows),
                "s_real": s_real,
                "s_imag": s_imag,
                "temporal_status": _temporal_status_for_root(
                    result.get("physical_root_metadata") or [],
                    s_real,
                    s_imag,
                ),
                "q_plus_real": q_plus_real,
                "q_plus_imag": q_plus_imag,
                "q_minus_real": q_minus_real,
                "q_minus_imag": q_minus_imag,
                "absolute_residual": verification.get("absolute_residual"),
                "relative_residual": verification.get("relative_residual"),
                "high_precision_refinement_status": (
                    (verification.get("high_precision_refinement") or {}).get("status")
                ),
                "admissibility_status": physical.get("status"),
            }
        )
    return rows


def _temporal_status_for_root(
    metadata_rows: Iterable[Mapping[str, Any]],
    s_real: float | None,
    s_imag: float | None,
) -> str | None:
    if s_real is None or s_imag is None:
        return None
    for row in metadata_rows:
        root = row.get("root")
        if root is None:
            continue
        real = float(root["real"])
        imag = float(root["imag"])
        if np.isclose(real, s_real, rtol=1e-10, atol=1e-12) and np.isclose(
            imag,
            s_imag,
            rtol=1e-10,
            atol=1e-12,
        ):
            return row.get("temporal_status")
    return None


def _point_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    points = [record for record in records if record.get("record_type") == "point"]
    indices = [int(record["index"]) for record in points]
    duplicates = sorted(index for index, count in Counter(indices).items() if count > 1)
    if duplicates:
        raise SweepSummaryError(
            "Duplicate point record(s) found: " + ", ".join(str(i) for i in duplicates)
        )
    return sorted(points, key=lambda row: int(row["index"]))


def _metadata_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    metadata = [record for record in records if record.get("record_type") == "metadata"]
    if not metadata:
        raise SweepSummaryError("No metadata record found")
    return metadata


def _infer_grid_parameter_order(definition: Mapping[str, Any]) -> tuple[str, ...]:
    swept_values = definition.get("swept_values") or {}
    names = tuple(swept_values)
    if len(names) <= 1:
        return names

    points = sorted(definition.get("points") or [], key=lambda row: int(row["index"]))
    if not points:
        return names

    values_by_name = {
        name: tuple(float(value) for value in swept_values[name])
        for name in names
    }
    for order in itertools.permutations(names):
        combinations = itertools.product(*(values_by_name[name] for name in order))
        expected = [
            {name: value for name, value in zip(order, combo, strict=True)}
            for combo in combinations
        ]
        if len(expected) != len(points):
            continue
        matched = True
        for point, combo in zip(points, expected, strict=True):
            coordinates = point.get("coordinates") or {}
            for name, value in combo.items():
                if not np.isclose(float(coordinates[name]), value, rtol=0.0, atol=1e-14):
                    matched = False
                    break
            if not matched:
                break
        if matched:
            return order
    return names


def _parameter_ordering(definition: Mapping[str, Any]) -> tuple[str, ...]:
    canonical = ("Arho", "Amu", "AG", "Lambda", "Ek")
    if definition.get("kind") == "cartesian_grid":
        fixed_values = definition.get("fixed_values") or {}
        fixed_names = tuple(name for name in canonical if name in fixed_values)
        return fixed_names + _infer_grid_parameter_order(definition)
    first_point = (definition.get("points") or [{}])[0]
    coordinates = first_point.get("coordinates") or {}
    return tuple(name for name in canonical if name in coordinates)


def read_sweep_summary(path: str | Path) -> SweepSummary:
    """Read a Sprint 6 JSONL file and build compact summary rows."""

    source = Path(path)
    records = _load_jsonl(source)
    metadata_records = _metadata_records(records)
    metadata = metadata_records[-1]
    definition = metadata.get("definition")
    if not isinstance(definition, dict):
        raise SweepSummaryError("Metadata record does not contain a sweep definition")

    point_records = _point_records(records)
    expected_count = int(definition.get("point_count", len(definition.get("points") or [])))
    if len(point_records) != expected_count:
        raise SweepSummaryError(
            f"Expected {expected_count} point record(s), found {len(point_records)}"
        )

    expected_indices = set(range(expected_count))
    actual_indices = {int(record["index"]) for record in point_records}
    if actual_indices != expected_indices:
        missing = sorted(expected_indices - actual_indices)
        extra = sorted(actual_indices - expected_indices)
        raise SweepSummaryError(f"Point index mismatch; missing={missing}, extra={extra}")

    point_rows = tuple(_extract_point_row(record) for record in point_records)
    root_rows: list[dict[str, Any]] = []
    for record in point_records:
        root_rows.extend(_extract_root_rows(record))

    validation = validate_summary_rows(point_rows, tuple(root_rows), expected_count)
    return SweepSummary(
        source_jsonl_path=str(source),
        workflow_version=metadata.get("workflow_version"),
        sweep_definition=definition,
        analysis_options=metadata.get("analyze_kwargs") or {},
        creation_time=metadata.get("created_unix_time"),
        parameter_ordering=_parameter_ordering(definition),
        coordinate_system=str(definition.get("coordinate_system", "unknown")),
        point_rows=point_rows,
        admissible_root_rows=tuple(root_rows),
        aggregate_classification_counts=_classification_counts(point_rows),
        validation_report=validation,
        metadata_records=tuple(metadata_records),
    )


def _classification_counts(point_rows: Iterable[Mapping[str, Any]]) -> dict[str, int]:
    counter = Counter(row["final_classification"] for row in point_rows)
    return {
        classification: int(counter.get(classification, 0))
        for classification in CLASSIFICATION_CODES
    }


def validate_summary_rows(
    point_rows: tuple[dict[str, Any], ...],
    root_rows: tuple[dict[str, Any], ...],
    expected_point_count: int,
) -> dict[str, Any]:
    """Validate internal consistency of extracted summary rows."""

    point_ids = [int(row["point_id"]) for row in point_rows]
    duplicate_ids = sorted(
        point_id for point_id, count in Counter(point_ids).items() if count > 1
    )
    if duplicate_ids:
        raise SweepSummaryError(
            "Duplicate summary point id(s): "
            + ", ".join(str(point_id) for point_id in duplicate_ids)
        )
    if len(point_rows) != expected_point_count:
        raise SweepSummaryError(
            f"Expected {expected_point_count} point row(s), found {len(point_rows)}"
        )

    root_count_by_point = Counter(int(row["point_id"]) for row in root_rows)
    mismatches = []
    for row in point_rows:
        point_id = int(row["point_id"])
        expected = int(row["physically_admissible_root_count"])
        actual = int(root_count_by_point.get(point_id, 0))
        if expected != actual:
            mismatches.append(
                {
                    "point_id": point_id,
                    "expected": expected,
                    "actual": actual,
                }
            )
    if mismatches:
        raise SweepSummaryError(f"Root-count mismatch: {mismatches}")

    return {
        "point_count_matches_definition": True,
        "duplicate_point_ids": [],
        "root_counts_match_point_rows": True,
        "classification_counts": _classification_counts(point_rows),
        "root_row_count": len(root_rows),
    }


def write_summary_outputs(summary: SweepSummary, output_dir: str | Path) -> dict[str, str]:
    """Write CSV/NPZ/JSON summary outputs and return created paths."""

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}

    point_csv = output / "point_summary.csv"
    _write_csv(point_csv, POINT_SUMMARY_FIELDS, summary.point_rows)
    paths["point_summary_csv"] = str(point_csv)

    root_csv = output / "admissible_roots.csv"
    _write_csv(root_csv, ROOT_SUMMARY_FIELDS, summary.admissible_root_rows)
    paths["admissible_roots_csv"] = str(root_csv)

    if summary.sweep_definition.get("kind") == "cartesian_grid":
        grid_path = output / "grid_data.npz"
        _write_grid_npz(summary, grid_path)
        paths["grid_data_npz"] = str(grid_path)

    metadata_path = output / "metadata.json"
    metadata_payload = summary_metadata(summary, paths)
    metadata_path.write_text(
        json.dumps(metadata_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    paths["metadata_json"] = str(metadata_path)

    report_path = output / "validation_report.txt"
    report_path.write_text(_format_validation_report(summary), encoding="utf-8")
    paths["validation_report"] = str(report_path)
    return paths


def _write_csv(
    path: Path,
    fieldnames: tuple[str, ...],
    rows: Iterable[Mapping[str, Any]],
) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in fieldnames})


def _write_grid_npz(summary: SweepSummary, path: Path) -> None:
    definition = summary.sweep_definition
    order = _infer_grid_parameter_order(definition)
    swept_values = definition.get("swept_values") or {}
    shape = tuple(len(swept_values[name]) for name in order)
    if not shape:
        return

    arrays: dict[str, Any] = {
        "grid_parameter_names": np.asarray(order),
        "classification_code": np.full(shape, -1, dtype=int),
        "effective_polynomial_degree": np.full(shape, -1, dtype=int),
        "polynomial_candidate_count": np.full(shape, -1, dtype=int),
        "domain_valid_candidate_count": np.full(shape, -1, dtype=int),
        "mathematically_genuine_root_count": np.full(shape, -1, dtype=int),
        "physically_admissible_root_count": np.full(shape, -1, dtype=int),
        "marginal_root_count": np.full(shape, -1, dtype=int),
        "physically_nonadmissible_root_count": np.full(shape, -1, dtype=int),
        "runtime_seconds": np.full(shape, np.nan, dtype=float),
    }
    for name in order:
        arrays[f"coord_{name}"] = np.asarray(swept_values[name], dtype=float)

    for row in summary.point_rows:
        multi_index = np.unravel_index(int(row["point_id"]), shape, order="C")
        arrays["classification_code"][multi_index] = int(row["classification_code"])
        degree = row["effective_polynomial_degree"]
        arrays["effective_polynomial_degree"][multi_index] = (
            -1 if degree in {None, ""} else int(degree)
        )
        for field in (
            "polynomial_candidate_count",
            "domain_valid_candidate_count",
            "mathematically_genuine_root_count",
            "physically_admissible_root_count",
            "marginal_root_count",
            "physically_nonadmissible_root_count",
        ):
            arrays[field][multi_index] = int(row[field])
        arrays["runtime_seconds"][multi_index] = float(row["runtime_seconds"])

    np.savez(path, **arrays)


def summary_metadata(summary: SweepSummary, output_paths: Mapping[str, str] | None = None) -> dict[str, Any]:
    """Return JSON-safe metadata for exported summary datasets."""

    return {
        "source_jsonl_path": summary.source_jsonl_path,
        "workflow_version": summary.workflow_version,
        "sweep_definition": summary.sweep_definition,
        "analysis_options": summary.analysis_options,
        "creation_time": summary.creation_time,
        "parameter_ordering": list(summary.parameter_ordering),
        "coordinate_system": summary.coordinate_system,
        "classification_codes": CLASSIFICATION_CODES,
        "classification_code_note": (
            "Categorical identifiers only; no ordinal meaning is implied."
        ),
        "aggregate_classification_counts": summary.aggregate_classification_counts,
        "validation_report": summary.validation_report,
        "output_paths": dict(output_paths or {}),
    }


def _format_validation_report(summary: SweepSummary) -> str:
    lines = [
        "Sweep summary validation report",
        "",
        f"Source JSONL: {summary.source_jsonl_path}",
        f"Workflow version: {summary.workflow_version}",
        f"Point rows: {len(summary.point_rows)}",
        f"Admissible root rows: {len(summary.admissible_root_rows)}",
        f"Coordinate system: {summary.coordinate_system}",
        f"Parameter ordering: {', '.join(summary.parameter_ordering)}",
        "",
        "Classification counts:",
    ]
    for classification, count in summary.aggregate_classification_counts.items():
        lines.append(f"  {classification}: {count}")
    lines.extend(["", "Validation:", json.dumps(summary.validation_report, indent=2)])
    return "\n".join(lines) + "\n"
