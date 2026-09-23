"""Compact, checksummed and atomic production result storage."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Iterable

import numpy as np


STORAGE_SCHEMA_VERSION = 1
STATUS_CODES = {
    "outside_theorem_scope": 0,
    "identically_zero_polynomial": 1,
    "constant_nonzero_polynomial": 2,
    "only_domain_invalid_candidates": 3,
    "only_spurious_candidates": 4,
    "genuine_roots_nonadmissible": 5,
    "marginal_spatial_decay": 6,
    "physically_admissible_roots": 7,
    "numerically_unresolved": 8,
}


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: str | Path, value: Any) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=target.name + ".", suffix=".tmp", dir=target.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, target)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def _root_rows(record: dict[str, Any], save_all_candidates: bool):
    result = record.get("result") or {}
    verification = result.get("verification_results") or []
    physical = result.get("physical_results") or []
    original = result.get("polynomial_candidates") or []
    rows = []
    for index, verification_row in enumerate(verification):
        physical_row = physical[index] if index < len(physical) else {}
        if not save_all_candidates and not physical_row.get("physically_admissible", False):
            continue
        final = verification_row.get("candidate") or {}
        initial = original[index] if index < len(original) else final
        rows.append((initial, final, verification_row, physical_row))
    return rows


def records_to_arrays(
    records: Iterable[dict[str, Any]], *, save_all_candidates: bool
) -> dict[str, np.ndarray]:
    rows = list(records)
    point_count = len(rows)
    offsets = [0]
    roots = []
    real_counts = []
    pair_counts = []
    unmatched_counts = []
    for record in rows:
        point_roots = _root_rows(record, save_all_candidates)
        roots.extend(point_roots)
        offsets.append(len(roots))
        admissible = [
            complex(float(root[1]["real"]), float(root[1]["imag"]))
            for root in point_roots if root[3].get("physically_admissible", False)
        ]
        pending = list(admissible)
        real_count = pair_count = unmatched_count = 0
        while pending:
            root = pending.pop(0)
            tolerance = 1e-10 + 1e-8 * max(1.0, abs(root))
            if abs(root.imag) <= tolerance:
                real_count += 1
                continue
            match = next((index for index, other in enumerate(pending)
                if abs(other - root.conjugate()) <= 1e-10 + 1e-8 * max(1.0, abs(root), abs(other))), None)
            if match is None:
                unmatched_count += 1
            else:
                pending.pop(match)
                pair_count += 1
        real_counts.append(real_count)
        pair_counts.append(pair_count)
        unmatched_counts.append(unmatched_count)

    def complex_parts(value):
        value = value or {}
        return float(value.get("real", np.nan)), float(value.get("imag", np.nan))

    def finite_or_nan(value):
        return np.nan if value is None else float(value)

    arrays: dict[str, np.ndarray] = {
        "schema_version": np.asarray(STORAGE_SCHEMA_VERSION, dtype=np.int16),
        "save_all_candidates": np.asarray(save_all_candidates, dtype=np.bool_),
        "point_index": np.asarray([row["index"] for row in rows], dtype=np.int64),
        "runtime_seconds": np.asarray([row.get("runtime_seconds", np.nan) for row in rows]),
        "classification": np.asarray(
            [STATUS_CODES[row["classification"]] for row in rows], dtype=np.uint8
        ),
        "error_flag": np.asarray([row.get("error") is not None for row in rows], dtype=np.bool_),
        "candidate_count": np.asarray([
            int((((row.get("result") or {}).get("counts") or {}).get("total_polynomial_candidates", 0)))
            for row in rows
        ], dtype=np.int16),
        "admissible_count": np.asarray([
            int((((row.get("result") or {}).get("counts") or {}).get("physically_admissible_roots", 0)))
            for row in rows
        ], dtype=np.int16),
        "admissible_real_count": np.asarray(real_counts, dtype=np.int16),
        "admissible_conjugate_pair_count": np.asarray(pair_counts, dtype=np.int16),
        "admissible_unmatched_nonreal_count": np.asarray(unmatched_counts, dtype=np.int16),
        "warning_count": np.asarray([
            len((row.get("result") or {}).get("warnings") or []) for row in rows
        ], dtype=np.int16),
        "high_precision_candidate_count": np.asarray([
            sum(
                verification.get("high_precision_refinement") is not None
                for verification in ((row.get("result") or {}).get("verification_results") or [])
            )
            for row in rows
        ], dtype=np.int16),
        "root_offsets": np.asarray(offsets, dtype=np.int64),
    }
    coordinate_names = ("Ek", "Arho", "Amu", "AG", "Lambda")
    for name in coordinate_names:
        arrays[f"coord_{name}"] = np.asarray([
            (row.get("coordinates") or {}).get(name, (row.get("parameters") or {}).get(name, np.nan))
            for row in rows
        ], dtype=np.float64)

    initial_parts = [complex_parts(root[0]) for root in roots]
    final_parts = [complex_parts(root[1]) for root in roots]
    arrays.update({
        "root_initial_real": np.asarray([x[0] for x in initial_parts]),
        "root_initial_imag": np.asarray([x[1] for x in initial_parts]),
        "root_real": np.asarray([x[0] for x in final_parts]),
        "root_imag": np.asarray([x[1] for x in final_parts]),
        "absolute_residual": np.asarray([finite_or_nan(x[2].get("absolute_residual")) for x in roots]),
        "relative_residual": np.asarray([finite_or_nan(x[2].get("relative_residual")) for x in roots]),
        "verification_threshold": np.asarray([finite_or_nan(x[2].get("verification_threshold")) for x in roots]),
        "domain_valid": np.asarray([x[2].get("domain_valid", False) for x in roots], dtype=np.bool_),
        "mathematically_genuine": np.asarray([x[2].get("mathematically_genuine", False) for x in roots], dtype=np.bool_),
        "physically_admissible": np.asarray([x[3].get("physically_admissible", False) for x in roots], dtype=np.bool_),
        "hp_precision": np.asarray([
            int((x[2].get("high_precision_refinement") or {}).get("working_precision", 0)) for x in roots
        ], dtype=np.int16),
        "verification_status": np.asarray([str(x[2].get("status") or "") for x in roots], dtype="U32"),
        "rejection_reason": np.asarray([str(x[2].get("rejection_reason") or "") for x in roots], dtype="U96"),
        "admissibility_status": np.asarray([str(x[3].get("status") or "") for x in roots], dtype="U32"),
    })
    for q_name in ("q_plus", "q_minus"):
        parts = [complex_parts(root[3].get(q_name)) for root in roots]
        arrays[f"{q_name}_real"] = np.asarray([x[0] for x in parts])
        arrays[f"{q_name}_imag"] = np.asarray([x[1] for x in parts])
    return arrays


def _validate_arrays(arrays: dict[str, np.ndarray], expected_points: int | None = None) -> None:
    if int(arrays["schema_version"]) != STORAGE_SCHEMA_VERSION:
        raise ValueError("unsupported storage schema")
    points = arrays["point_index"]
    offsets = arrays["root_offsets"]
    if expected_points is not None and len(points) != expected_points:
        raise ValueError("stored point count does not match expected count")
    if len(offsets) != len(points) + 1 or offsets[0] != 0:
        raise ValueError("invalid root offsets")
    if np.any(np.diff(offsets) < 0) or offsets[-1] != len(arrays["root_real"]):
        raise ValueError("root offsets do not match root arrays")
    if len(set(int(value) for value in points)) != len(points):
        raise ValueError("duplicate point indices")


def write_npz_atomic(path: str | Path, arrays: dict[str, np.ndarray], *, compressed: bool = True) -> str:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    _validate_arrays(arrays)
    fd, temporary_name = tempfile.mkstemp(prefix=target.name + ".", suffix=".npz", dir=target.parent)
    os.close(fd)
    try:
        if compressed:
            np.savez_compressed(temporary_name, **arrays)
        else:
            np.savez(temporary_name, **arrays)
        with np.load(temporary_name, allow_pickle=False) as loaded:
            _validate_arrays({name: loaded[name] for name in loaded.files})
        os.replace(temporary_name, target)
        return sha256_file(target)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def write_summary_atomic(path: str | Path, payload_paths: Iterable[str | Path]) -> None:
    point_index = []
    classification = []
    candidate_count = []
    admissible_count = []
    admissible_real_count = []
    admissible_conjugate_pair_count = []
    admissible_unmatched_nonreal_count = []
    runtime_seconds = []
    high_precision_candidate_count = []
    for payload_path in payload_paths:
        with np.load(payload_path, allow_pickle=False) as loaded:
            point_index.append(loaded["point_index"])
            classification.append(loaded["classification"])
            candidate_count.append(loaded["candidate_count"])
            admissible_count.append(loaded["admissible_count"])
            admissible_real_count.append(loaded["admissible_real_count"])
            admissible_conjugate_pair_count.append(loaded["admissible_conjugate_pair_count"])
            admissible_unmatched_nonreal_count.append(loaded["admissible_unmatched_nonreal_count"])
            runtime_seconds.append(loaded["runtime_seconds"])
            high_precision_candidate_count.append(loaded["high_precision_candidate_count"])
    order = np.argsort(np.concatenate(point_index))
    arrays = {
        "schema_version": np.asarray(STORAGE_SCHEMA_VERSION, dtype=np.int16),
        "point_index": np.concatenate(point_index)[order],
        "classification": np.concatenate(classification)[order],
        "candidate_count": np.concatenate(candidate_count)[order],
        "admissible_count": np.concatenate(admissible_count)[order],
        "admissible_real_count": np.concatenate(admissible_real_count)[order],
        "admissible_conjugate_pair_count": np.concatenate(admissible_conjugate_pair_count)[order],
        "admissible_unmatched_nonreal_count": np.concatenate(admissible_unmatched_nonreal_count)[order],
        "runtime_seconds": np.concatenate(runtime_seconds)[order],
        "high_precision_candidate_count": np.concatenate(high_precision_candidate_count)[order],
    }
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary_name = tempfile.mkstemp(prefix=target.name + ".", suffix=".npz", dir=target.parent)
    os.close(fd)
    try:
        np.savez_compressed(temporary_name, **arrays)
        with np.load(temporary_name, allow_pickle=False) as loaded:
            if len(loaded["point_index"]) != len(set(loaded["point_index"].tolist())):
                raise ValueError("summary contains duplicate point indices")
        os.replace(temporary_name, target)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def validate_payload(path: str | Path, sidecar: str | Path, *, config_sha256: str, case_id: str) -> bool:
    data_path = Path(path)
    sidecar_path = Path(sidecar)
    if not data_path.is_file() or not sidecar_path.is_file():
        return False
    try:
        metadata = json.loads(sidecar_path.read_text(encoding="utf-8"))
        if metadata.get("schema_version") != STORAGE_SCHEMA_VERSION:
            return False
        if metadata.get("complete") is not True:
            return False
        # Failed points must be recomputed when a case is selected from a retry
        # manifest; otherwise a checksum-valid, incomplete scientific result
        # would be skipped forever.
        if int(metadata.get("failure_count", 0)) != 0:
            return False
        if metadata.get("config_sha256") != config_sha256 or metadata.get("case_id") != case_id:
            return False
        if metadata.get("sha256") != sha256_file(data_path):
            return False
        with np.load(data_path, allow_pickle=False) as loaded:
            arrays = {name: loaded[name] for name in loaded.files}
        _validate_arrays(arrays, int(metadata["point_count"]))
        return True
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return False


def write_payload(
    path: str | Path,
    sidecar: str | Path,
    records: Iterable[dict[str, Any]],
    *,
    case_id: str,
    config_sha256: str,
    save_all_candidates: bool,
    bounds: list[int] | None = None,
    compressed: bool = True,
) -> dict[str, Any]:
    rows = list(records)
    arrays = records_to_arrays(rows, save_all_candidates=save_all_candidates)
    checksum = write_npz_atomic(path, arrays, compressed=compressed)
    metadata = {
        "schema_version": STORAGE_SCHEMA_VERSION,
        "case_id": case_id,
        "config_sha256": config_sha256,
        "point_count": len(rows),
        "bounds": bounds,
        "sha256": checksum,
        "complete": True,
        "failure_count": sum(row.get("error") is not None for row in rows),
    }
    atomic_json(sidecar, metadata)
    return metadata


def tile_bounds(shape: tuple[int, int], tile_shape: tuple[int, int]):
    rows, columns = shape
    tile_rows, tile_columns = tile_shape
    for row_start in range(0, rows, tile_rows):
        for column_start in range(0, columns, tile_columns):
            yield (
                row_start,
                min(rows, row_start + tile_rows),
                column_start,
                min(columns, column_start + tile_columns),
            )
