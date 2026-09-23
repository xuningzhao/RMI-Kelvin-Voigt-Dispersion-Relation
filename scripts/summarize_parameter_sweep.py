#!/usr/bin/env python
"""Create compact summary datasets from a Sprint 6 sweep JSONL file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.sweep_summary import read_sweep_summary, write_summary_outputs  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_jsonl", type=Path, help="Sprint 6 JSONL sweep output")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for point/root/grid/metadata summary outputs",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir
    if output_dir is None:
        output_dir = args.input_jsonl.with_suffix("").parent / (
            args.input_jsonl.with_suffix("").name + "_summary"
        )

    summary = read_sweep_summary(args.input_jsonl)
    paths = write_summary_outputs(summary, output_dir)

    print("Sweep summary")
    print(f"  source       = {summary.source_jsonl_path}")
    print(f"  points       = {len(summary.point_rows)}")
    print(f"  roots        = {len(summary.admissible_root_rows)}")
    print(f"  coordinate   = {summary.coordinate_system}")
    print(f"  output dir   = {output_dir}")
    print()
    print("Classification counts")
    for name, count in summary.aggregate_classification_counts.items():
        print(f"  {name:32s} {count}")
    print()
    print("Created files")
    for name, path in paths.items():
        print(f"  {name:24s} {path}")


if __name__ == "__main__":
    main()

