"""Global production result index assembled from case-local status files."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .production_storage import atomic_json


def build_global_index(campaign_root: str | Path) -> dict[str, Any]:
    root = Path(campaign_root)
    records = []
    for path in sorted(root.glob("group_*/**/status.json")):
        try:
            record = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            record = {
                "case_id": path.parent.name,
                "status": "corrupt_status",
                "complete": False,
                "failure_information": f"could not read {path}",
                "output_path": str(path.parent),
            }
        records.append(record)
    records.sort(key=lambda row: (str(row.get("group", "")), str(row.get("case_id", ""))))
    index = {
        "schema_version": 1,
        "case_count": len(records),
        "complete_case_count": sum(bool(row.get("complete")) for row in records),
        "failed_case_count": sum(int(row.get("failure_count", 0)) > 0 for row in records),
        "cases": records,
    }
    atomic_json(root / "global_index.json", index)
    return index
