#!/usr/bin/env python
"""Map the primary complex-to-real Ek transition over (Arho, Amu) at fixed AG."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.run_mp_reference_sweep import _point as _mp_point  # noqa: E402
from src.parameter_sweep import lambda_from_ek  # noqa: E402
from src.workflow import analyze_parameter_point  # noqa: E402


def _down_crossings(states: np.ndarray) -> np.ndarray:
    """Indices whose interval changes from complex to non-complex."""
    values = np.asarray(states, dtype=bool)
    return np.flatnonzero(values[:-1] & ~values[1:])


def _float_complex(arho: float, amu: float, ag: float, ek: float, imag_tol: float) -> bool:
    result = analyze_parameter_point(
        arho, amu, ag, lambda_from_ek(ek), enable_high_precision_refinement=False
    )
    roots = (
        result.physical_run.physically_admissible_roots
        if result.physical_run is not None
        else ()
    )
    return any(abs(complex(root).imag) > imag_tol for root in roots)


def _mp_complex(
    arho: float, amu: float, ag: float, ek: float, precision: int, imag_tol: float
) -> bool:
    payload = (0, repr(arho), repr(amu), repr(ag), repr(ek), precision)
    accepted = _mp_point(payload)[3]
    return any(abs(root[0].imag) > imag_tol for root in accepted)


def _locate(payload: tuple) -> dict:
    (
        point_id, arho, amu, ag, ek_min, ek_max, scan_samples,
        precision, imag_tol, ek_tolerance,
    ) = payload
    ek = np.linspace(ek_min, ek_max, scan_samples)
    coarse = np.asarray(
        [_float_complex(arho, amu, ag, float(value), imag_tol) for value in ek],
        dtype=bool,
    )
    coarse_crossings = _down_crossings(coarse)
    if not len(coarse_crossings):
        status = "always_complex" if np.all(coarse) else (
            "always_noncomplex" if not np.any(coarse) else "no_down_crossing"
        )
        return {
            "point_id": point_id, "Arho": arho, "Amu": amu,
            "Ek_c": None, "status": status, "coarse_crossing_count": 0,
            "mp_iterations": 0, "bracket_width": None,
        }

    # The primary transition is the last complex-to-real crossing as Ek increases.
    # Search candidate brackets from right to left in case float64 and MP disagree.
    selected = None
    mp_evaluations = 0
    for crossing in reversed(coarse_crossings):
        left, right = float(ek[crossing]), float(ek[crossing + 1])
        left_state = _mp_complex(arho, amu, ag, left, precision, imag_tol)
        right_state = _mp_complex(arho, amu, ag, right, precision, imag_tol)
        mp_evaluations += 2
        if left_state and not right_state:
            selected = (left, right)
            break
    if selected is None:
        # Backend disagreement can shift or create a crossing.  Fall back to
        # the complete coarse MP scan for this point rather than dropping it.
        mp_coarse = np.asarray(
            [
                _mp_complex(arho, amu, ag, float(value), precision, imag_tol)
                for value in ek
            ],
            dtype=bool,
        )
        mp_evaluations += len(ek)
        mp_crossings = _down_crossings(mp_coarse)
        if not len(mp_crossings):
            return {
                "point_id": point_id, "Arho": arho, "Amu": amu,
                "Ek_c": None, "status": "mp_no_down_crossing",
                "coarse_crossing_count": int(len(coarse_crossings)),
                "mp_iterations": mp_evaluations, "bracket_width": None,
            }
        crossing = int(mp_crossings[-1])
        selected = (float(ek[crossing]), float(ek[crossing + 1]))

    left, right = selected
    while right - left > ek_tolerance:
        middle = 0.5 * (left + right)
        if _mp_complex(arho, amu, ag, middle, precision, imag_tol):
            left = middle
        else:
            right = middle
        mp_evaluations += 1
    return {
        "point_id": point_id, "Arho": arho, "Amu": amu,
        "Ek_c": 0.5 * (left + right), "status": "transition",
        "coarse_crossing_count": int(len(coarse_crossings)),
        "mp_iterations": mp_evaluations, "bracket_width": right - left,
    }


def _completed(path: Path) -> dict[int, dict]:
    rows: dict[int, dict] = {}
    if not path.exists():
        return rows
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            try:
                row = json.loads(line)
                rows[int(row["point_id"])] = row
            except (json.JSONDecodeError, KeyError, TypeError, ValueError):
                continue
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=51)
    parser.add_argument("--contrast-limit", type=float, default=0.999998)
    parser.add_argument("--AG", type=float, default=-1.0)
    parser.add_argument("--ek-min", type=float, default=-0.99)
    parser.add_argument("--ek-max", type=float, default=0.99)
    parser.add_argument("--scan-samples", type=int, default=65)
    parser.add_argument("--ek-tolerance", type=float, default=2.5e-4)
    parser.add_argument("--imag-tolerance", type=float, default=2.0e-3)
    parser.add_argument("--precision", type=int, default=80)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--progress-interval", type=int, default=100)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint = args.output_dir / "transition_points.jsonl"
    if checkpoint.exists() and not args.resume:
        parser.error(f"checkpoint exists: {checkpoint}; pass --resume or choose a new output")
    rows = _completed(checkpoint) if args.resume else {}
    contrasts = np.linspace(-args.contrast_limit, args.contrast_limit, args.samples)
    shape = (args.samples, args.samples)
    payloads = (
        (
            i * args.samples + j, float(arho), float(amu), args.AG,
            args.ek_min, args.ek_max, args.scan_samples, args.precision,
            args.imag_tolerance, args.ek_tolerance,
        )
        for i, arho in enumerate(contrasts)
        for j, amu in enumerate(contrasts)
        if i * args.samples + j not in rows
    )
    with checkpoint.open("a", encoding="utf-8") as handle:
        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            for completed, row in enumerate(
                executor.map(_locate, payloads, chunksize=1), start=1
            ):
                handle.write(json.dumps(row, separators=(",", ":")) + "\n")
                handle.flush()
                rows[int(row["point_id"])] = row
                if completed % args.progress_interval == 0:
                    print(f"completed_this_run={completed}; total={len(rows)}/{np.prod(shape)}", flush=True)

    ekc = np.full(shape, np.nan)
    status = np.full(shape, "missing", dtype="U24")
    crossing_count = np.zeros(shape, dtype=np.int8)
    for point_id, row in rows.items():
        index = np.unravel_index(point_id, shape)
        if row["Ek_c"] is not None:
            ekc[index] = float(row["Ek_c"])
        status[index] = row["status"]
        crossing_count[index] = int(row["coarse_crossing_count"])
    np.savez_compressed(
        args.output_dir / "transition_map.npz",
        coord_Arho=contrasts, coord_Amu=contrasts, Ek_critical=ekc,
        transition_status=status, coarse_crossing_count=crossing_count,
    )
    metadata = {
        "fixed_AG": args.AG, "shape": list(shape),
        "contrast_range": [-args.contrast_limit, args.contrast_limit],
        "Ek_scan_range": [args.ek_min, args.ek_max],
        "coarse_scan_samples": args.scan_samples,
        "coarse_backend": "float64", "refinement_backend": "mp",
        "mp_decimal_precision": args.precision,
        "imaginary_tolerance": args.imag_tolerance,
        "Ek_target_bracket_width": args.ek_tolerance,
        "primary_transition": "last complex-to-real crossing with increasing Ek",
    }
    (args.output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )
    finite = np.isfinite(ekc)
    print(f"transitions={np.count_nonzero(finite)}/{ekc.size}")
    if np.any(finite):
        print(f"Ek,c range=[{np.nanmin(ekc):.8f}, {np.nanmax(ekc):.8f}]")


if __name__ == "__main__":
    main()
