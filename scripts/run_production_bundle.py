#!/usr/bin/env python
"""Run one production job with process-based case workers; never submits Slurm."""

from __future__ import annotations

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.production_config import load_production_config  # noqa: E402
from src.production_manifest import read_manifest  # noqa: E402
from src.production_runner import run_case  # noqa: E402
from src.production_storage import atomic_json  # noqa: E402


def _run_lane(payload):
    records, config_path, campaign_root, maximum_new_units = payload
    config = load_production_config(config_path)
    statuses = []
    for case in records:
        try:
            statuses.append(run_case(case, config, campaign_root,
                maximum_new_units=maximum_new_units))
        except Exception as exc:  # noqa: BLE001 - isolate case failures
            root = Path(campaign_root) / case["output_path"]
            status = {
                "case_id": case["case_id"], "group": case["group"],
                "status": "case_worker_failed", "complete": False,
                "failure_count": 1, "error": type(exc).__name__,
                "error_message": str(exc), "output_path": str(root),
                "config_sha256": config.sha256,
                "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
            }
            atomic_json(root / "status.json", status)
            statuses.append(status)
    return statuses


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--campaign-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--case-workers", type=int, required=True)
    parser.add_argument("--cases-per-worker", type=int, required=True)
    parser.add_argument("--maximum-new-units", type=int)
    args = parser.parse_args()
    config = load_production_config(args.config)
    allocated_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
    if allocated_cpus is not None and args.case_workers > int(allocated_cpus):
        raise SystemExit("--case-workers exceeds SLURM_CPUS_PER_TASK")
    records = read_manifest(args.manifest)
    expected = args.case_workers * args.cases_per_worker
    if len(records) != expected:
        raise SystemExit(f"manifest has {len(records)} cases; expected {expected}")
    lanes = tuple(records[index:index + args.cases_per_worker]
        for index in range(0, len(records), args.cases_per_worker))
    payloads = tuple((lane, args.config, args.campaign_root, args.maximum_new_units)
        for lane in lanes)
    if args.case_workers == 1:
        nested = [_run_lane(payloads[0])]
    else:
        with ProcessPoolExecutor(max_workers=args.case_workers) as executor:
            nested = list(executor.map(_run_lane, payloads))
    statuses = [status for lane in nested for status in lane]
    print(json.dumps(statuses, indent=2, sort_keys=True))
    if any(not status.get("complete") or int(status.get("failure_count", 0)) > 0
        for status in statuses):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
