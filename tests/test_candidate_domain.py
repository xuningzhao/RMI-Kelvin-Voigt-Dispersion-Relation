"""Tests for polynomial-candidate definition-domain checks."""

from __future__ import annotations

import unittest

import numpy as np

from src.candidate_domain import (
    candidate_intermediates,
    check_candidate_domain,
    check_candidates_domain,
)
from src.dispersion import D_star_Lambda
from src.polynomial import candidate_roots, effective_polynomial_from_parameters


class CandidateDomainTests(unittest.TestCase):
    def test_generic_complex_candidate_is_domain_valid(self):
        result = check_candidate_domain(0.37 + 0.91j, 0.2, -0.3, 0.4, 1.7)

        self.assertTrue(result.domain_valid)
        self.assertEqual(result.failed_conditions, ())
        self.assertIsNotNone(result.intermediates)
        self.assertTrue(all(check.passed for check in result.checks))

    def test_polynomial_candidates_receive_domain_diagnostics(self):
        effective = effective_polynomial_from_parameters(0.2, -0.3, 0.4, 1.7)
        roots = candidate_roots(effective)
        results = check_candidates_domain(roots, 0.2, -0.3, 0.4, 1.7)

        self.assertEqual(len(results), roots.size)
        self.assertTrue(any(result.domain_valid for result in results))
        for result in results:
            self.assertTrue(result.checks)

    def test_s_near_zero_fails_s_nonzero_condition(self):
        result = check_candidate_domain(1e-16 + 0j, 0.0, 0.0, 0.0, 1.0)

        self.assertFalse(result.domain_valid)
        self.assertIn("s", result.failed_conditions)

    def test_E_plus_zero_fails_E_plus_condition(self):
        Arho, Amu, AG, Lambda = 0.1, 0.2, 0.1, 1.0
        s = -Lambda * (1.0 + AG) / (1.0 + Amu)

        result = check_candidate_domain(s, Arho, Amu, AG, Lambda)

        self.assertFalse(result.domain_valid)
        self.assertIn("E_plus", result.failed_conditions)
        self.assertIn("skipped_radicals_after_E_denominator_failure", result.notes)

    def test_E_minus_zero_fails_E_minus_condition(self):
        Arho, Amu, AG, Lambda = 0.1, 0.2, 0.1, 1.0
        s = -Lambda * (1.0 - AG) / (1.0 - Amu)

        result = check_candidate_domain(s, Arho, Amu, AG, Lambda)

        self.assertFalse(result.domain_valid)
        self.assertIn("E_minus", result.failed_conditions)
        self.assertIn("skipped_radicals_after_E_denominator_failure", result.notes)

    def test_B_minus_zero_fails_B_minus_condition(self):
        # With Arho=1, q_minus=sqrt(1)=1.  At s=-Lambda and Amu=0,
        # E_plus=-E_minus, so B_minus=E_plus+E_minus*q_minus=0.
        result = check_candidate_domain(-1.0, 1.0, 0.0, 0.4, 1.0)

        self.assertFalse(result.domain_valid)
        self.assertIn("B_minus", result.failed_conditions)
        self.assertIsNotNone(result.intermediates)
        self.assertAlmostEqual(abs(result.intermediates.B_minus), 0.0)

    def test_B_plus_zero_fails_B_plus_condition(self):
        # With Arho=-1, q_plus=sqrt(1)=1.  At s=-Lambda and Amu=0,
        # E_minus=-E_plus, so B_plus=E_minus+E_plus*q_plus=0.
        result = check_candidate_domain(-1.0, -1.0, 0.0, 0.4, 1.0)

        self.assertFalse(result.domain_valid)
        self.assertIn("B_plus", result.failed_conditions)
        self.assertIsNotNone(result.intermediates)
        self.assertAlmostEqual(abs(result.intermediates.B_plus), 0.0)

    def test_principal_numpy_square_root_is_used(self):
        intermediates = candidate_intermediates(-0.4 + 1.3j, -0.2, 0.35, -0.5, 2.1)

        self.assertEqual(intermediates.q_plus, complex(np.sqrt(intermediates.R_plus)))
        self.assertEqual(intermediates.q_minus, complex(np.sqrt(intermediates.R_minus)))

    def test_intermediates_are_consistent_with_dispersion_denominators(self):
        s = -0.4 + 0.8j
        Arho, Amu, AG, Lambda = 0.25, -0.35, 0.6, 1.4
        intermediates = candidate_intermediates(s, Arho, Amu, AG, Lambda)

        eta1 = (1.0 + Amu) + (Lambda / s) * (1.0 + AG)
        eta2 = (1.0 - Amu) + (Lambda / s) * (1.0 - AG)
        q1 = np.sqrt(1.0 + ((1.0 + Arho) * s) / eta1)
        q2 = np.sqrt(1.0 + ((1.0 - Arho) * s) / eta2)
        denominator1 = eta1 + eta2 * q2
        denominator2 = eta2 + eta1 * q1

        self.assertAlmostEqual(intermediates.E_plus / s, eta1)
        self.assertAlmostEqual(intermediates.E_minus / s, eta2)
        self.assertAlmostEqual(intermediates.q_plus, q1)
        self.assertAlmostEqual(intermediates.q_minus, q2)
        self.assertAlmostEqual(intermediates.B_minus / s, denominator1)
        self.assertAlmostEqual(intermediates.B_plus / s, denominator2)

        # If the domain check passes, src.dispersion must be evaluable at the
        # same point.  This is not a root verification; it only checks
        # definition consistency.
        result = check_candidate_domain(s, Arho, Amu, AG, Lambda)
        self.assertTrue(result.domain_valid)
        value = D_star_Lambda(s, Arho, Amu, AG, Lambda)
        self.assertTrue(np.isfinite(value.real))
        self.assertTrue(np.isfinite(value.imag))


if __name__ == "__main__":
    unittest.main()

