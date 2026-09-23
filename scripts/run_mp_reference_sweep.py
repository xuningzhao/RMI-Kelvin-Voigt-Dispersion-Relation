#!/usr/bin/env python
"""Arbitrary-precision reference sweep using the unchanged canonical P14."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import os
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.mp_reference import solve_reference_point  # noqa: E402


def _point(payload):
    index, arho_text, amu_text, ag_text, ek_text, precision = payload
    result = solve_reference_point(
        arho_text, amu_text, ag_text, ek_text, decimal_precision=precision,
    )
    return (index, result.polynomial_candidate_count,
            result.mathematically_genuine_count, list(result.accepted))


def _modal(real_count, total):
    result = np.zeros(total.shape, dtype=np.int8)
    has_real = real_count > 0
    has_complex = total > real_count
    result[has_real & ~has_complex] = 1
    result[~has_real & has_complex] = 2
    result[has_real & has_complex] = 3
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=501)
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--precision", type=int, default=80)
    parser.add_argument("--Arho", default="0")
    parser.add_argument("--Amu", default="0")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--progress-interval", type=int, default=5000)
    args = parser.parse_args()
    # Preserve the command-line decimal strings all the way into mpmath.
    arho_text = str(args.Arho)
    amu_text = str(args.Amu)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ag = np.linspace(-1.0, 1.0, args.samples)
    ek = np.linspace(-0.99, 0.99, args.samples)
    shape = (len(ag), len(ek))
    count = np.zeros(shape, dtype=np.int16)
    real_count = np.zeros(shape, dtype=np.int16)
    genuine_count = np.zeros(shape, dtype=np.int16)
    candidate_count = np.zeros(shape, dtype=np.int16)
    rows = []
    payloads = (
        (i * len(ek) + j, arho_text, amu_text,
         repr(float(a)), repr(float(e)), args.precision)
        for i, a in enumerate(ag) for j, e in enumerate(ek)
    )
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for completed, (index, candidates, genuine, accepted) in enumerate(
            executor.map(_point, payloads, chunksize=4), start=1
        ):
            i, j = np.unravel_index(index, shape)
            candidate_count[i, j] = candidates
            genuine_count[i, j] = genuine
            count[i, j] = len(accepted)
            real_count[i, j] = sum(abs(root[0].imag) <= 1e-10 for root in accepted)
            for root_index, (root, qp, qm, residual) in enumerate(accepted):
                rows.append((index, root_index, ag[i], ek[j], root.real, root.imag,
                             qp.real, qp.imag, qm.real, qm.imag, residual))
            if completed % args.progress_interval == 0:
                print(f"completed={completed}/{count.size}", flush=True)
    modal = _modal(real_count, count)
    np.savez_compressed(
        args.output_dir / "grid_data.npz", coord_AG=ag, coord_Ek=ek,
        physically_admissible_root_count=count,
        admissible_real_root_count=real_count,
        modal_configuration=modal,
        mathematically_genuine_root_count=genuine_count,
        polynomial_candidate_count=candidate_count,
    )
    with (args.output_dir / "admissible_roots.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(("point_id", "root_index", "AG", "Ek", "s_real", "s_imag",
                         "q_plus_real", "q_plus_imag", "q_minus_real", "q_minus_imag",
                         "high_precision_relative_residual"))
        writer.writerows(rows)
    metadata = {
        "backend": "mp", "decimal_precision": args.precision,
        "coefficient_source": "src/polynomial_coefficients.py canonical expressions",
        "lambda_evaluation": "mp: ((1-Ek)/(1+Ek))**2",
        "physical_admissibility": {
            "criterion": "Re(q) > atol + rtol * max(1, abs(q))",
            "q_relative_tolerance": 1e-10,
            "q_absolute_tolerance": 1e-12,
        },
        "fixed_parameters": {"Arho": arho_text, "Amu": amu_text},
        "shape": list(shape), "workers": args.workers,
    }
    (args.output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
