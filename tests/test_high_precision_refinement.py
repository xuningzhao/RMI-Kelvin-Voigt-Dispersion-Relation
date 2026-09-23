"""Tests for targeted high-precision refinement of ambiguous candidates."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import mpmath as mp

from src.high_precision_refinement import (
    mp_all_polynomial_roots,
    mp_lambda_from_ek,
    mp_p14_coefficients_from_ek,
)
from src.parameter_sweep import lambda_from_ek
from src.polynomial_backend import polynomial_roots_from_ek
from src.workflow import analyze_parameter_point


ARHO = -0.95
AMU = -0.8


class HighPrecisionRefinementTests(unittest.TestCase):
    def test_mp_backend_computes_lambda_directly_from_ek(self):
        value = mp_lambda_from_ek("0.90", decimal_precision=80)
        self.assertLess(abs(value - mp.mpf(1) / 361), mp.mpf("1e-79"))

    def test_mp_backend_retains_large_ek_slow_root_candidates(self):
        coefficients, Lambda = mp_p14_coefficients_from_ek(
            "0", "0", "-1", "0.90", decimal_precision=80
        )
        roots = mp_all_polynomial_roots(
            coefficients, decimal_precision=80, max_steps=2000
        )
        nearby = sorted(roots, key=lambda root: abs(root + Lambda))[:2]
        self.assertEqual(len(roots), 12)
        self.assertLess(abs(nearby[0] + Lambda), mp.mpf("2e-6"))

    def test_selectable_mp_backend_uses_mp_lambda(self):
        roots, value = polynomial_roots_from_ek(
            "0", "0", "-1", "0.90", backend="mp", decimal_precision=80
        )
        self.assertEqual(len(roots), 12)
        self.assertIsInstance(value, mp.mpf)

    def test_mp_backend_deflates_exact_zero_roots_with_multiplicity(self):
        coefficients, _ = mp_p14_coefficients_from_ek(
            "-0.05", "0", "0", "-0.99", decimal_precision=80
        )
        roots = mp_all_polynomial_roots(coefficients, decimal_precision=80)
        self.assertEqual(len(roots), 14)
        self.assertEqual(sum(root == 0 for root in roots), 8)

    def test_previously_unresolved_point_recovers_extra_genuine_roots(self):
        Lambda = lambda_from_ek(0.6933333333333334)

        baseline = analyze_parameter_point(ARHO, AMU, 0.7466666666666668, Lambda)
        refined = analyze_parameter_point(
            ARHO,
            AMU,
            0.7466666666666668,
            Lambda,
            enable_high_precision_refinement=True,
        )

        self.assertEqual(baseline.classification, "numerically_unresolved")
        self.assertGreater(baseline.counts.ambiguous_mathematical_candidates, 0)
        self.assertEqual(refined.classification, "physically_admissible_roots")
        self.assertEqual(refined.counts.ambiguous_mathematical_candidates, 0)
        self.assertGreater(
            refined.counts.mathematically_genuine_roots,
            baseline.counts.mathematically_genuine_roots,
        )
        self.assertGreater(
            refined.counts.physically_admissible_roots,
            baseline.counts.physically_admissible_roots,
        )

    def test_clustered_two_root_case_preserves_refined_roots(self):
        refined = analyze_parameter_point(
            ARHO,
            AMU,
            0.7466666666666668,
            lambda_from_ek(0.6933333333333334),
            enable_high_precision_refinement=True,
        )
        hp_rows = [
            row.high_precision_refinement
            for row in refined.verification_run.results
            if row.high_precision_refinement is not None
        ]

        self.assertEqual(len(hp_rows), 2)
        self.assertTrue(all(row.status == "refined_genuine" for row in hp_rows))
        self.assertTrue(all(row.cluster_size >= 2 for row in hp_rows))
        refined_roots = [row.refined_candidate for row in hp_rows]
        self.assertEqual(len(set(refined_roots)), len(refined_roots))

    def test_robust_spurious_point_does_not_invoke_high_precision(self):
        result = analyze_parameter_point(
            ARHO,
            AMU,
            -0.32,
            lambda_from_ek(-0.8),
            enable_high_precision_refinement=True,
        )

        self.assertEqual(result.classification, "only_spurious_candidates")
        self.assertEqual(result.counts.mathematically_genuine_roots, 0)
        self.assertFalse(
            any(
                row.high_precision_refinement is not None
                for row in result.verification_run.results
            )
        )

    def test_normal_point_does_not_invoke_high_precision(self):
        result = analyze_parameter_point(
            0.2,
            -0.3,
            0.4,
            1.7,
            enable_high_precision_refinement=True,
        )

        self.assertEqual(result.classification, "physically_admissible_roots")
        self.assertFalse(
            any(
                row.high_precision_refinement is not None
                for row in result.verification_run.results
            )
        )

    def test_refinement_failure_remains_unresolved(self):
        with patch(
            "src.high_precision_refinement.mp_p14_coefficients",
            side_effect=RuntimeError("forced refinement failure"),
        ):
            result = analyze_parameter_point(
                ARHO,
                AMU,
                0.7466666666666668,
                lambda_from_ek(0.6933333333333334),
                enable_high_precision_refinement=True,
            )

        self.assertEqual(result.classification, "numerically_unresolved")
        self.assertGreater(result.counts.ambiguous_mathematical_candidates, 0)
        hp_rows = [
            row.high_precision_refinement
            for row in result.verification_run.results
            if row.high_precision_refinement is not None
        ]
        self.assertTrue(hp_rows)
        self.assertTrue(all(row.status == "refinement_failed" for row in hp_rows))

    def test_refined_roots_keep_conjugate_structure_for_complex_case(self):
        result = analyze_parameter_point(
            ARHO,
            AMU,
            -0.8,
            lambda_from_ek(-0.8),
            enable_high_precision_refinement=True,
        )
        roots = [
            row.candidate
            for row in result.verification_run.results
            if row.mathematically_genuine
        ]

        self.assertEqual(len(roots), 2)
        for root in roots:
            self.assertTrue(any(abs(other - root.conjugate()) < 1e-8 for other in roots))


if __name__ == "__main__":
    unittest.main()
