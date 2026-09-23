"""Tests for mathematical verification of P14 candidates."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np

from src.candidate_domain import check_candidate_domain
from src.dispersion import D_star_Lambda
from src.polynomial import candidate_roots, effective_polynomial_from_parameters
from src.root_verification import (
    residual_scale_from_domain,
    verification_threshold,
    verify_candidate_root,
    verify_candidate_roots,
    verify_p14_roots_from_parameters,
)


GENERIC_PARAMS = (0.2, -0.3, 0.4, 1.7)
REDUCED_PARAMS = (0.3, -0.3, 0.2, 0.8)
NONSYMMETRIC_PARAMS = (-0.65, 0.15, -0.5, 2.25)


class RootVerificationTests(unittest.TestCase):
    def test_known_genuine_roots_are_identified(self):
        run = verify_p14_roots_from_parameters(*GENERIC_PARAMS)

        self.assertEqual(run.summary.total_polynomial_candidates, 14)
        self.assertEqual(run.summary.domain_valid_candidates, 14)
        self.assertEqual(run.summary.mathematically_genuine_roots, 4)
        self.assertEqual(len(run.genuine_roots), 4)

        for result in run.results:
            if result.mathematically_genuine:
                self.assertLessEqual(
                    result.absolute_residual,
                    result.verification_threshold,
                )
                self.assertIsNone(result.rejection_reason)

    def test_perturbed_genuine_root_fails_verification(self):
        run = verify_p14_roots_from_parameters(*GENERIC_PARAMS)
        root = run.genuine_roots[0]
        perturbed = root + 1e-3 * (1.0 + 0.5j)

        result = verify_candidate_root(perturbed, *GENERIC_PARAMS)

        self.assertFalse(result.mathematically_genuine)
        self.assertIn(result.status, {"residual_too_large", "numerically_ambiguous"})
        self.assertIsNotNone(result.absolute_residual)

    def test_domain_invalid_candidate_is_not_evaluated(self):
        with patch("src.root_verification.D_star_Lambda") as evaluator:
            result = verify_candidate_root(1e-16, 0.0, 0.0, 0.0, 1.0)

        evaluator.assert_not_called()
        self.assertFalse(result.domain_valid)
        self.assertFalse(result.mathematically_genuine)
        self.assertEqual(result.status, "domain_invalid")
        self.assertIsNone(result.absolute_residual)
        self.assertIn("s", result.domain_result.failed_conditions)

    def test_reduced_degree_case_verifies_candidates(self):
        effective = effective_polynomial_from_parameters(*REDUCED_PARAMS)
        self.assertEqual(effective.status, "reduced_degree")
        self.assertEqual(effective.effective_degree, 13)

        run = verify_p14_roots_from_parameters(*REDUCED_PARAMS)
        self.assertEqual(run.summary.total_polynomial_candidates, 13)
        self.assertEqual(run.summary.mathematically_genuine_roots, 4)

    def test_complex_conjugate_genuine_roots_are_present_for_real_parameters(self):
        run = verify_p14_roots_from_parameters(*NONSYMMETRIC_PARAMS)
        roots = np.asarray(run.genuine_roots)

        self.assertEqual(roots.size, 4)
        for root in roots:
            distances = np.abs(roots - np.conj(root))
            self.assertLess(np.min(distances), 1e-8)

    def test_residual_values_match_D_star_Lambda_evaluator(self):
        effective = effective_polynomial_from_parameters(*GENERIC_PARAMS)
        roots = candidate_roots(effective)
        run = verify_candidate_roots(roots, *GENERIC_PARAMS)

        for result in run.results:
            if not result.domain_valid:
                continue
            expected = abs(D_star_Lambda(result.candidate, *GENERIC_PARAMS))
            self.assertAlmostEqual(result.absolute_residual, expected, places=14)

    def test_scale_aware_threshold_matches_documented_policy(self):
        candidate = -2.540599781523809 + 1.1850371393599597j
        domain = check_candidate_domain(candidate, *GENERIC_PARAMS)
        scale = residual_scale_from_domain(domain)
        threshold = verification_threshold(
            scale,
            relative_tolerance=1e-8,
            absolute_tolerance=1e-10,
        )

        result = verify_candidate_root(candidate, *GENERIC_PARAMS, domain_result=domain)

        self.assertAlmostEqual(result.residual_scale, scale)
        self.assertAlmostEqual(result.verification_threshold, threshold)
        self.assertAlmostEqual(
            result.relative_residual,
            result.absolute_residual / max(1.0, scale),
        )


if __name__ == "__main__":
    unittest.main()

