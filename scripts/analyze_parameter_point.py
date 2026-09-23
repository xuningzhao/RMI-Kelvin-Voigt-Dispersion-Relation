#!/usr/bin/env python
"""Analyze one nondimensional material-parameter tuple.

This script runs the finalized single-point workflow only.  It does not perform
parameter sweeps, plotting, or phase-diagram generation.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.workflow import analyze_parameter_point  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Arho", type=float, required=True)
    parser.add_argument("--Amu", type=float, required=True)
    parser.add_argument("--AG", type=float, required=True)
    parser.add_argument("--Lambda", type=float, required=True)
    parser.add_argument(
        "--diagnostic-outside-theorem-scope",
        action="store_true",
        help="Run the numerical pipeline outside theorem scope as diagnostic only.",
    )
    parser.add_argument(
        "--enable-high-precision-refinement",
        action="store_true",
        help="Refine only numerically ambiguous candidates using high precision.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    result = analyze_parameter_point(
        args.Arho,
        args.Amu,
        args.AG,
        args.Lambda,
        analyze_outside_theorem_scope=args.diagnostic_outside_theorem_scope,
        enable_high_precision_refinement=args.enable_high_precision_refinement,
    )

    print("Single parameter-point analysis")
    print()
    print("Input parameters")
    print(f"  Arho   = {result.parameters.Arho:.17g}")
    print(f"  Amu    = {result.parameters.Amu:.17g}")
    print(f"  AG     = {result.parameters.AG:.17g}")
    print(f"  Lambda = {result.parameters.Lambda:.17g}")
    print()
    print("Theorem scope")
    print(f"  applies           = {result.theorem_scope.applies}")
    print(f"  failed_conditions = {','.join(result.theorem_scope.failed_conditions) or 'none'}")
    print(f"  analysis_skipped  = {result.theorem_scope.analysis_skipped}")
    print(f"  diagnostic_only   = {result.theorem_scope.diagnostic_only}")
    print()
    print("Polynomial metadata")
    if result.effective_polynomial is None:
        print("  not evaluated")
    else:
        print(f"  status           = {result.effective_polynomial.status}")
        print(f"  effective_degree = {result.effective_polynomial.effective_degree}")
        print(f"  reduced_degree   = {result.reduced_degree}")
        print(f"  coefficient_scale = {result.effective_polynomial.coefficient_scale:.6e}")
    print()
    print("Counts")
    counts = result.counts
    print(f"  polynomial_candidates          = {counts.total_polynomial_candidates}")
    print(f"  domain_valid_candidates        = {counts.domain_valid_candidates}")
    print(f"  domain_invalid_candidates      = {counts.domain_invalid_candidates}")
    print(f"  mathematically_genuine_roots   = {counts.mathematically_genuine_roots}")
    print(f"  ambiguous_math_candidates      = {counts.ambiguous_mathematical_candidates}")
    print(f"  physically_admissible_roots    = {counts.physically_admissible_roots}")
    print(f"  marginal_spatial_decay_roots   = {counts.marginal_spatial_decay_roots}")
    print(f"  physically_nonadmissible_roots = {counts.physically_nonadmissible_roots}")
    print(f"  evaluation_failures            = {counts.evaluation_failures}")
    print()
    print(f"Final classification: {result.classification}")
    if result.warnings:
        print("Warnings")
        for warning in result.warnings:
            print(f"  - {warning}")
    print()

    if result.verification_run is not None:
        print("Mathematically genuine roots")
        genuine = [
            row for row in result.verification_run.results if row.mathematically_genuine
        ]
        if not genuine:
            print("  none")
        for index, row in enumerate(genuine):
            print(
                f"  genuine[{index:02d}] = {row.candidate.real:.17e}"
                f" {row.candidate.imag:+.17e}j"
                f"    |D*| = {row.absolute_residual:.6e}"
                f"    rel_D = {row.relative_residual:.6e}"
            )

    if result.physical_run is not None:
        print()
        print("Physically admissible roots")
        accepted = [
            row for row in result.physical_run.results if row.physically_admissible
        ]
        if not accepted:
            print("  none")
        for index, row in enumerate(accepted):
            print(
                f"  physical[{index:02d}] = {row.candidate.real:.17e}"
                f" {row.candidate.imag:+.17e}j"
                f"    Re(q+) = {row.re_q_plus:.6e}"
                f"    Re(q-) = {row.re_q_minus:.6e}"
                f"    reason = {row.reason}"
            )

    if result.physical_root_metadata:
        print()
        print("Temporal metadata for physically admissible roots")
        for index, row in enumerate(result.physical_root_metadata):
            print(
                f"  root[{index:02d}] temporal_status = {row.temporal_status}"
                f"    Re(s) = {row.real_part:.6e}"
                f"    Im(s) = {row.imaginary_part:.6e}"
            )


if __name__ == "__main__":
    main()
