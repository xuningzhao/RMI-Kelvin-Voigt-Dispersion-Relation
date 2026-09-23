#!/usr/bin/env python
"""Build deterministic Group A and Group B production manifests."""

from __future__ import annotations

import argparse
import json
import hashlib
import platform
from pathlib import Path
import subprocess
import sys

import numpy as np
import mpmath

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.production_config import load_production_config  # noqa: E402
from src.production_manifest import GROUPS, build_cases, write_manifest  # noqa: E402
from src.production_storage import atomic_json  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "config/production_sweep_v1.json")
    parser.add_argument("--campaign-root", type=Path, required=True)
    args = parser.parse_args()
    config = load_production_config(args.config)
    manifest_dir = args.campaign_root / "manifests"
    counts = {}
    for group in GROUPS:
        cases = build_cases(config, group)
        name = group.lower() + ".jsonl"
        write_manifest(manifest_dir / name, cases)
        counts[group] = len(cases)
    def command(*parts: str) -> str | None:
        try:
            return subprocess.run(parts, cwd=ROOT, check=True, text=True,
                capture_output=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return None

    diff = command("git", "diff", "--no-ext-diff", "--binary")
    provenance = {
        "git_commit": command("git", "rev-parse", "HEAD"),
        "git_status": command("git", "status", "--short"),
        "git_diff_sha256": None if not diff else hashlib.sha256(diff.encode()).hexdigest(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "mpmath_version": mpmath.__version__,
        "platform": platform.platform(),
    }
    atomic_json(args.campaign_root / "configuration.json", {
        "configuration": config.data, "config_sha256": config.sha256,
        "source": str(args.config), "manifest_counts": counts,
        "provenance": provenance,
    })
    print(json.dumps(counts, sort_keys=True))


if __name__ == "__main__":
    main()
