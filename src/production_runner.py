"""Bundled production execution with Group A cases and Group B tiles."""

from __future__ import annotations

import json
import os
from pathlib import Path
import time
import traceback
from typing import Any

from .parameter_sweep import SweepPointInput, SweepPointResult, sweep_point_to_json
from .production_config import ProductionConfig
from .production_storage import (
    atomic_json,
    tile_bounds,
    validate_payload,
    write_payload,
    write_summary_atomic,
)
from .workflow import analyze_parameter_point


def _analyze(payload: tuple[SweepPointInput, dict[str, Any]]) -> dict[str, Any]:
    point, options = payload
    start = time.perf_counter()
    try:
        result = analyze_parameter_point(
            point.parameters["Arho"], point.parameters["Amu"],
            point.parameters["AG"], point.parameters["Lambda"], **options
        )
        row = SweepPointResult(point.index, point.coordinates, point.parameters,
            result, result.classification, time.perf_counter() - start)
    except Exception as exc:  # noqa: BLE001 - point isolation is required
        row = SweepPointResult(point.index, point.coordinates, point.parameters,
            None, "numerically_unresolved", time.perf_counter() - start,
            error=type(exc).__name__, traceback_summary="".join(
                traceback.format_exception_only(type(exc), exc)).strip())
    return sweep_point_to_json(row)


def _serial_map(case: dict[str, Any], config: ProductionConfig, bounds=None):
    options = _solver_options(config)
    for point in _case_points(case, bounds):
        yield _analyze((point, options))


def _point(index: int, coordinates: dict[str, float]) -> SweepPointInput:
    from .parameter_sweep import lambda_from_ek
    coordinates = dict(coordinates)
    lam = lambda_from_ek(coordinates["Ek"])
    coordinates["Lambda"] = lam
    parameters = {name: float(coordinates[name]) for name in ("Arho", "Amu", "AG")}
    parameters["Lambda"] = lam
    return SweepPointInput(index, coordinates, parameters)


def _case_points(case: dict[str, Any], bounds: tuple[int, int, int, int] | None = None):
    swept = tuple(case["swept"])
    axes = case["axes"]
    fixed = {key: float(value) for key, value in case["fixed"].items()}
    if len(swept) == 1:
        for index, value in enumerate(axes[swept[0]]):
            yield _point(index, {**fixed, swept[0]: float(value)})
        return
    shape = tuple(len(axes[name]) for name in swept)
    row_start, row_stop, column_start, column_stop = bounds or (0, shape[0], 0, shape[1])
    for row in range(row_start, row_stop):
        for column in range(column_start, column_stop):
            index = row * shape[1] + column
            yield _point(index, {**fixed, swept[0]: float(axes[swept[0]][row]), swept[1]: float(axes[swept[1]][column])})


def _solver_options(config: ProductionConfig) -> dict[str, Any]:
    return config.solver_options


def run_case(
    case: dict[str, Any],
    config: ProductionConfig,
    campaign_root: str | Path,
    *,
    maximum_new_units: int | None = None,
) -> dict[str, Any]:
    if case.get("config_sha256") != config.sha256:
        raise ValueError("case manifest and production configuration hashes differ")
    root = Path(campaign_root) / case["output_path"]
    root.mkdir(parents=True, exist_ok=True)
    if int(config.data["execution"]["point_workers"]) != 1:
        raise ValueError("production cases must run serially within case workers")
    save_all = config.save_all_candidates
    compressed = bool(config.data["storage"]["compression"])
    start = time.perf_counter()
    checkpoint_seconds = 0.0
    failures = 0
    failure_rows: list[dict[str, Any]] = []
    payload_paths: list[Path] = []
    completed_units = 0
    skipped_units = 0
    atomic_json(root / "metadata.json", {"case": case, "configuration": config.data})

    if case["group"] == "A":
        data_path, sidecar = root / "spectrum.npz", root / "spectrum.json"
        if validate_payload(data_path, sidecar, config_sha256=config.sha256, case_id=case["case_id"]):
            skipped_units = 1
            failures += int(json.loads(sidecar.read_text(encoding="utf-8"))["failure_count"])
        else:
            rows = list(_serial_map(case, config))
            failure_rows.extend({"index": row["index"], "error": row["error"], "traceback_summary": row["traceback_summary"]} for row in rows if row.get("error"))
            checkpoint_start = time.perf_counter()
            metadata = write_payload(data_path, sidecar, rows, case_id=case["case_id"],
                config_sha256=config.sha256, save_all_candidates=save_all, compressed=compressed)
            checkpoint_seconds += time.perf_counter() - checkpoint_start
            failures += metadata["failure_count"]
            completed_units = 1
        payload_paths.append(data_path)
        expected_units = 1
    else:
        shape = tuple(int(x) for x in case["shape"])
        tile_shape = tuple(int(x) for x in config.data["execution"]["group_b_tile_shape"])
        bounds_list = list(tile_bounds(shape, tile_shape))
        expected_units = len(bounds_list)
        tile_dir = root / "tiles"
        for tile_index, bounds in enumerate(bounds_list):
            stem = f"tile_{tile_index:04d}_r{bounds[0]:03d}-{bounds[1]:03d}_c{bounds[2]:03d}-{bounds[3]:03d}"
            data_path, sidecar = tile_dir / f"{stem}.npz", tile_dir / f"{stem}.json"
            if validate_payload(data_path, sidecar, config_sha256=config.sha256, case_id=case["case_id"]):
                skipped_units += 1
                failures += int(json.loads(sidecar.read_text(encoding="utf-8"))["failure_count"])
                payload_paths.append(data_path)
                continue
            if maximum_new_units is not None and completed_units >= maximum_new_units:
                break
            rows = list(_serial_map(case, config, bounds))
            failure_rows.extend({"index": row["index"], "error": row["error"], "traceback_summary": row["traceback_summary"]} for row in rows if row.get("error"))
            checkpoint_start = time.perf_counter()
            metadata = write_payload(data_path, sidecar, rows, case_id=case["case_id"],
                config_sha256=config.sha256, save_all_candidates=save_all,
                bounds=list(bounds), compressed=compressed)
            checkpoint_seconds += time.perf_counter() - checkpoint_start
            failures += metadata["failure_count"]
            completed_units += 1
            payload_paths.append(data_path)

    summary_start = time.perf_counter()
    write_summary_atomic(root / "summary.npz", payload_paths)
    checkpoint_seconds += time.perf_counter() - summary_start
    import numpy as np
    with np.load(root / "summary.npz", allow_pickle=False) as summary:
        stored_points = len(summary["point_index"])
        point_work_seconds = float(summary["runtime_seconds"].sum())
        high_precision_candidates = int(summary["high_precision_candidate_count"].sum())
        unresolved_points = int((summary["classification"] == 8).sum())
    complete = stored_points == int(case["point_count"])
    if maximum_new_units is None and not complete:
        raise RuntimeError("completed payloads do not cover the expected case points")
    if failure_rows:
        atomic_json(root / "failures.json", {"case_id": case["case_id"], "failures": failure_rows})
    else:
        (root / "failures.json").unlink(missing_ok=True)

    status = {
        "case_id": case["case_id"], "group": case["group"],
        "status": ("incomplete" if not complete else
            ("complete" if failures == 0 else "complete_with_failures")),
        "complete": complete, "runtime_seconds": time.perf_counter() - start,
        "expected_units": expected_units, "completed_units": completed_units,
        "skipped_units": skipped_units, "failure_count": failures,
        "output_path": str(root), "config_sha256": config.sha256,
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "updated_unix_time": time.time(),
        "stored_points": stored_points,
        "point_work_seconds": point_work_seconds,
        "checkpoint_seconds": checkpoint_seconds,
        "high_precision_candidates": high_precision_candidates,
        "unresolved_points": unresolved_points,
        "output_bytes": sum(path.stat().st_size for path in root.rglob("*") if path.is_file()),
    }
    atomic_json(root / "status.json", status)
    return status
