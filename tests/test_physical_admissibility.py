"""Tests for physical admissibility of mathematically genuine roots."""

from __future__ import annotations

import unittest

from src.candidate_domain import (
    CandidateDomainResult,
    CandidateIntermediates,
)
from src.physical_admissibility import (
    assess_physical_admissibility,
    assess_physical_admissibility_many,
    classify_q_decay,
    q_real_tolerance,
)
from src.root_verification import (
    CandidateVerificationResult,
    verify_p14_roots_from_parameters,
)


GENERIC_PARAMS = (0.2, -0.3, 0.4, 1.7)


def _synthetic_verification(q_plus: complex, q_minus: complex) -> CandidateVerificationResult:
    intermediates = CandidateIntermediates(
        s=1.0 + 0.25j,
        E_plus=2.0 + 0.1j,
        E_minus=1.5 - 0.2j,
        R_plus=q_plus**2,
        R_minus=q_minus**2,
        q_plus=q_plus,
        q_minus=q_minus,
        B_plus=1.0 + 0.1j,
        B_minus=1.2 - 0.1j,
    )
    domain = CandidateDomainResult(
        s=intermediates.s,
        domain_valid=True,
        failed_conditions=(),
        intermediates=intermediates,
        checks=(),
        notes=(),
    )
    return CandidateVerificationResult(
        candidate=intermediates.s,
        domain_result=domain,
        domain_valid=True,
        absolute_residual=0.0,
        relative_residual=0.0,
        residual_scale=1.0,
        verification_threshold=1e-10,
        mathematically_genuine=True,
        status="mathematically_genuine",
        rejection_reason=None,
        dispersion_value=0.0,
    )


class PhysicalAdmissibilityTests(unittest.TestCase):
    def test_physically_admissible_genuine_roots_are_identified(self):
        verification = verify_p14_roots_from_parameters(*GENERIC_PARAMS)
        physical = assess_physical_admissibility_many(verification.results)

        self.assertEqual(physical.summary.total_candidates, 14)
        self.assertEqual(physical.summary.mathematically_genuine_roots, 4)
        self.assertEqual(physical.summary.physically_admissible_roots, 4)
        self.assertEqual(len(physical.physically_admissible_roots), 4)

    def test_q_values_are_reused_from_domain_diagnostics(self):
        verification = verify_p14_roots_from_parameters(*GENERIC_PARAMS)
        genuine = next(result for result in verification.results if result.mathematically_genuine)
        physical = assess_physical_admissibility(genuine)

        self.assertIsNotNone(genuine.domain_result.intermediates)
        self.assertEqual(physical.q_plus, genuine.domain_result.intermediates.q_plus)
        self.assertEqual(physical.q_minus, genuine.domain_result.intermediates.q_minus)
        self.assertEqual(physical.re_q_plus, physical.q_plus.real)
        self.assertEqual(physical.re_q_minus, physical.q_minus.real)

    def test_failed_upper_decay_is_classified(self):
        status, reason, accepted, _tol_plus, _tol_minus = classify_q_decay(
            -1.0 + 0.2j,
            1.0 + 0.1j,
        )

        self.assertEqual(status, "physically_nonadmissible")
        self.assertEqual(reason, "failed_upper_decay")
        self.assertFalse(accepted)

    def test_failed_lower_decay_is_classified(self):
        status, reason, accepted, _tol_plus, _tol_minus = classify_q_decay(
            1.0 + 0.2j,
            -1.0 + 0.1j,
        )

        self.assertEqual(status, "physically_nonadmissible")
        self.assertEqual(reason, "failed_lower_decay")
        self.assertFalse(accepted)

    def test_failed_both_decay_is_classified(self):
        status, reason, accepted, _tol_plus, _tol_minus = classify_q_decay(
            -1.0 + 0.2j,
            -1.0 + 0.1j,
        )

        self.assertEqual(status, "physically_nonadmissible")
        self.assertEqual(reason, "failed_both_decay")
        self.assertFalse(accepted)

    def test_near_spatial_decay_boundary_is_marginal(self):
        q_plus = 1e-14 + 2.0j
        q_minus = 1.0 + 0.5j
        result = assess_physical_admissibility(
            _synthetic_verification(q_plus, q_minus),
            absolute_tolerance=1e-12,
            relative_tolerance=0.0,
        )

        self.assertFalse(result.physically_admissible)
        self.assertEqual(result.status, "marginal_spatial_decay")
        self.assertEqual(result.reason, "marginal_spatial_decay")

    def test_complex_root_result_reports_complex_q_values(self):
        verification = verify_p14_roots_from_parameters(*GENERIC_PARAMS)
        genuine = next(
            result
            for result in verification.results
            if result.mathematically_genuine and abs(result.candidate.imag) > 0.1
        )
        physical = assess_physical_admissibility(genuine)

        self.assertIsInstance(physical.q_plus, complex)
        self.assertIsInstance(physical.q_minus, complex)
        self.assertGreater(physical.re_q_plus, physical.tolerance_plus)
        self.assertGreater(physical.re_q_minus, physical.tolerance_minus)

    def test_not_mathematically_genuine_roots_are_skipped(self):
        verification = verify_p14_roots_from_parameters(*GENERIC_PARAMS)
        rejected = next(result for result in verification.results if not result.mathematically_genuine)
        physical = assess_physical_admissibility(rejected)

        self.assertFalse(physical.physically_admissible)
        self.assertEqual(physical.status, "not_mathematically_genuine")
        self.assertIsNone(physical.q_plus)
        self.assertIsNone(physical.q_minus)

    def test_q_real_tolerance_is_scale_aware(self):
        small = q_real_tolerance(0.5 + 0.1j, relative_tolerance=1e-10, absolute_tolerance=1e-12)
        large = q_real_tolerance(10.0 + 0.1j, relative_tolerance=1e-10, absolute_tolerance=1e-12)

        self.assertGreater(large, small)


if __name__ == "__main__":
    unittest.main()

