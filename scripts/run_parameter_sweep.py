#!/usr/bin/env python
"""Run a small batch parameter sweep.

Examples
--------

Sweep AG and Lambda at fixed Arho/Amu:

    python scripts/run_parameter_sweep.py \\
      --Arho 0.2 --Amu -0.3 \\
      --sweep AG -0.5 0.0 0.5 \\
      --sweep Lambda 0.8 1.7

Sweep AG and Ek, recording both Ek and converted Lambda:

    python scripts/run_parameter_sweep.py \\
      --Arho 0.2 --Amu -0.3 \\
      --sweep AG -0.5 0.0 0.5 \\
      --sweep Ek -0.2 0.2 0.8
"""

from __future__ import annotations

import argparse
from datetime import datetime
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.parameter_sweep import cartesian_grid_definition, run_parameter_sweep  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Arho", type=float)
    parser.add_argument("--Amu", type=float)
    parser.add_argument("--AG", type=float)
    parser.add_argument("--Lambda", type=float)
    parser.add_argument("--Ek", type=float)
    parser.add_argument(
        "--sweep",
        action="append",
        nargs="+",
        metavar=("NAME", "VALUE"),
        help="Parameter sweep values, e.g. --sweep AG -0.5 0 0.5",
    )
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--diagnostic-outside-theorem-scope",
        action="store_true",
        help="Run outside-scope points diagnostically instead of skipping.",
    )
    parser.add_argument(
        "--enable-high-precision-refinement",
        action="store_true",
        help="Refine only numerically ambiguous candidates using high precision.",
    )
    parser.add_argument("--hp-initial-precision", type=int, default=80)
    parser.add_argument("--hp-maximum-precision", type=int, default=120)
    parser.add_argument("--hp-precision-step", type=int, default=20)
    return parser.parse_args()


def _parse_sweeps(raw_sweeps: list[list[str]] | None) -> dict[str, list[float]]:
    sweep: dict[str, list[float]] = {}
    for entry in raw_sweeps or []:
        if len(entry) < 2:
            raise SystemExit("--sweep requires a parameter name and at least one value")
        name = entry[0]
        if name not in {"Arho", "Amu", "AG", "Lambda", "Ek"}:
            raise SystemExit(f"Unsupported sweep parameter: {name}")
        if name in sweep:
            raise SystemExit(f"Duplicate sweep parameter: {name}")
        sweep[name] = [float(value) for value in entry[1:]]
    return sweep


def main() -> None:
    args = parse_args()
    sweep = _parse_sweeps(args.sweep)

    fixed = {}
    for name in ("Arho", "Amu", "AG", "Lambda", "Ek"):
        value = getattr(args, name)
        if value is not None and name not in sweep:
            fixed[name] = value

    if "Lambda" in fixed and "Ek" in fixed:
        raise SystemExit("Provide only one of fixed Lambda or fixed Ek")
    if "Lambda" in sweep and "Ek" in sweep:
        raise SystemExit("Sweep only one of Lambda or Ek")

    definition = cartesian_grid_definition(fixed=fixed, sweep=sweep)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = args.output or PROJECT_ROOT / "results" / f"parameter_sweep_{timestamp}.jsonl"

    print("Parameter sweep")
    print(f"  points  = {len(definition.points)}")
    print(f"  workers = {args.workers}")
    print(f"  output  = {output}")
    print()

    def progress(point, counts):
        done = sum(counts.values())
        unresolved = counts.get("numerically_unresolved", 0)
        print(
            f"[{done:04d}/{len(definition.points):04d}] "
            f"index={point.index} classification={point.classification} "
            f"unresolved={unresolved}"
        )

    result = run_parameter_sweep(
        definition,
        max_workers=args.workers,
        checkpoint_path=output,
        resume=args.resume,
        progress_callback=progress,
        analyze_outside_theorem_scope=args.diagnostic_outside_theorem_scope,
        enable_high_precision_refinement=args.enable_high_precision_refinement,
        high_precision_initial_decimal_precision=args.hp_initial_precision,
        high_precision_maximum_decimal_precision=args.hp_maximum_precision,
        high_precision_precision_step=args.hp_precision_step,
    )

    print()
    print("Final classification counts")
    for name, count in result.aggregate_classification_counts.items():
        print(f"  {name:32s} {count}")
    print()
    print(f"Unresolved points: {len(result.unresolved_points)}")
    print(f"Total runtime: {result.total_runtime_seconds:.3f} s")
    print(f"Output location: {result.checkpoint_path}")


if __name__ == "__main__":
    main()
