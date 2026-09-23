#!/usr/bin/env python
"""Validate production payload checksums and expected case completion."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.production_config import load_production_config  # noqa: E402
from src.production_manifest import read_manifest  # noqa: E402
from src.production_storage import tile_bounds, validate_payload  # noqa: E402


def validate_case(case, config, campaign_root):
    root = campaign_root / case["output_path"]
    invalid = []
    if case["group"] == "A":
        expected = [(root / "spectrum.npz", root / "spectrum.json")]
    else:
        shape = tuple(case["shape"])
        tile_shape = tuple(config.data["execution"]["group_b_tile_shape"])
        expected = []
        for index, bounds in enumerate(tile_bounds(shape, tile_shape)):
            stem = f"tile_{index:04d}_r{bounds[0]:03d}-{bounds[1]:03d}_c{bounds[2]:03d}-{bounds[3]:03d}"
            expected.append((root / "tiles" / f"{stem}.npz", root / "tiles" / f"{stem}.json"))
    for data, sidecar in expected:
        if not validate_payload(data, sidecar, config_sha256=config.sha256, case_id=case["case_id"]):
            invalid.append(str(data))
    return {"case_id": case["case_id"], "valid": not invalid, "invalid_payloads": invalid}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--campaign-root", type=Path, required=True)
    args = parser.parse_args()
    config = load_production_config(args.config)
    results = [validate_case(case, config, args.campaign_root) for case in read_manifest(args.manifest)]
    report = {"valid": all(row["valid"] for row in results), "cases": results}
    print(json.dumps(report, indent=2, sort_keys=True))
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
