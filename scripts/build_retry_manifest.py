#!/usr/bin/env python
"""Select incomplete or failed cases into a deterministic retry manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--manifest", type=Path, required=True)
parser.add_argument("--campaign-root", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
selected = []
with args.manifest.open(encoding="utf-8") as handle:
    for line in handle:
        case = json.loads(line)
        status_path = args.campaign_root / case["output_path"] / "status.json"
        try:
            status = json.loads(status_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            status = {}
        if not status.get("complete") or int(status.get("failure_count", 0)) > 0:
            selected.append(case)
args.output.parent.mkdir(parents=True, exist_ok=True)
temporary = args.output.with_suffix(args.output.suffix + ".tmp")
temporary.write_text("".join(json.dumps(case, sort_keys=True) + "\n" for case in selected), encoding="utf-8")
temporary.replace(args.output)
print(json.dumps({"retry_case_count": len(selected), "output": str(args.output)}))
