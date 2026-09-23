#!/usr/bin/env python3
"""Print a verified arbitrary-precision reference solve as JSON."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.mp_reference import solve_reference_point


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("Arho", "Amu", "AG", "Ek"):
        parser.add_argument(f"--{name}", required=True)
    parser.add_argument("--precision", type=int, default=80)
    args = parser.parse_args()
    parameters = {name: getattr(args, name) for name in ("Arho", "Amu", "AG", "Ek")}
    result = solve_reference_point(**parameters, decimal_precision=args.precision)
    def pair(value):
        return {"real": value.real, "imag": value.imag}
    print(json.dumps({
        "parameters": parameters,
        "backend": "mp",
        "decimal_precision": args.precision,
        "output_precision": "binary64; acceptance evaluated before conversion",
        "polynomial_candidate_count": result.polynomial_candidate_count,
        "mathematically_genuine_count": result.mathematically_genuine_count,
        "physically_admissible_count": len(result.accepted),
        "accepted": [
            {"s": pair(s), "q_plus": pair(qp), "q_minus": pair(qm),
             "relative_source_residual": residual}
            for s, qp, qm, residual in result.accepted
        ],
    }, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
