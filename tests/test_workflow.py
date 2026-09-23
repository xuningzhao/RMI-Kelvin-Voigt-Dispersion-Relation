"""Tests for single-parameter-point workflow orchestration."""

from __future__ import annotations

import unittest
from unittest.mock import patch

import numpy as np

from src.parameter_sweep import lambda_from_ek
from src.physical_admissibility import assess_physical_admissibility_many
from src.polynomial import build_effective_polynomial, candidate_roots
from src.root_verification import verify_candidate_roots
from src.workflow import (
    TheoremScopeStatus,
    WorkflowCounts,
    analyze_parameter_point,
    classify_parameter_point,
)


GENERIC_PARAMS = (0.2, -0.3, 0.4, 1.7)
ONLY_SPURIOUS_PARAMS = (0.5, 0.5, 0.0, 81.0)
REDUCED_PARAMS = (0.3, -0.3, 0.2, 0.8)
NEGATIVE_EK_REGRESSION_CASES = (
    (-0.8, -0.8, -0.44209312324515126 + 8.451250773470967j),
    (-0.64, -0.64, -0.446661767207 + 4.13308125598j),
    (0.64, -0.8, -2.74832226154 - 11.3650958193j),
    (0.64, -0.64, -1.62000870959 + 5.30724165312j),
)


def _scope() -> TheoremScopeStatus:
    return TheoremScopeStatus(
        applies=True,
        failed_conditions=(),
        analysis_skipped=False,
        diagnostic_only=False,
    )


def _counts(**overrides) -> WorkflowCounts:
    values = dict(
        effective_polynomial_degree=14,
        total_polynomial_candidates=14,
        domain_valid_candidates=14,
        domain_invalid_candidates=0,
        mathematically_genuine_roots=0,
        ambiguous_mathematical_candidates=0,
        physically_admissible_roots=0,
        marginal_spatial_decay_roots=0,
        physically_nonadmissible_roots=0,
        evaluation_failures=0,
    )
    values.update(overrides)
    return WorkflowCounts(**values)


class WorkflowTests(unittest.TestCase):
    def test_end_to_end_physically_admissible_classification(self):
        result = analyze_parameter_point(*GENERIC_PARAMS)

        self.assertTrue(result.theorem_scope.applies)
        self.assertEqual(result.classification, "physically_admissible_roots")
        self.assertEqual(result.counts.total_polynomial_candidates, 14)
        self.assertEqual(result.counts.mathematically_genuine_roots, 4)
        self.assertEqual(result.counts.physically_admissible_roots, 4)
        self.assertEqual(len(result.physical_root_metadata), 4)

    def test_end_to_end_only_spurious_classification(self):
        result = analyze_parameter_point(*ONLY_SPURIOUS_PARAMS)

        self.assertEqual(result.classification, "only_spurious_candidates")
        self.assertGreater(result.counts.domain_valid_candidates, 0)
        self.assertEqual(result.counts.mathematically_genuine_roots, 0)
        self.assertEqual(result.counts.physically_admissible_roots, 0)

    def test_reduced_degree_is_metadata_not_failure(self):
        result = analyze_parameter_point(*REDUCED_PARAMS)

        self.assertTrue(result.reduced_degree)
        self.assertEqual(result.counts.effective_polynomial_degree, 13)
        self.assertEqual(result.classification, "physically_admissible_roots")

    def test_negative_Ek_high_lambda_regressions_recover_physical_roots(self):
        for AG, Ek, reference_root in NEGATIVE_EK_REGRESSION_CASES:
            with self.subTest(AG=AG, Ek=Ek):
                result = analyze_parameter_point(-0.95, -0.8, AG, lambda_from_ek(Ek))

                self.assertEqual(result.effective_polynomial.status, "generic_degree")
                self.assertEqual(result.counts.effective_polynomial_degree, 14)
                self.assertEqual(result.counts.total_polynomial_candidates, 14)
                self.assertGreater(result.counts.mathematically_genuine_roots, 0)
                self.assertGreater(result.counts.physically_admissible_roots, 0)
                self.assertEqual(result.classification, "physically_admissible_roots")

                nearest_index = min(
                    range(len(result.polynomial_candidates)),
                    key=lambda index: abs(result.polynomial_candidates[index] - reference_root),
                )
                self.assertLess(
                    abs(result.polynomial_candidates[nearest_index] - reference_root),
                    1e-8,
                )
                self.assertTrue(result.domain_results[nearest_index].domain_valid)
                self.assertTrue(
                    result.verification_run.results[nearest_index].mathematically_genuine
                )
                self.assertTrue(
                    result.physical_run.results[nearest_index].physically_admissible
                )

    def test_outside_theorem_scope_is_skipped_by_default(self):
        result = analyze_parameter_point(0.2, -0.3, 0.4, 0.0)

        self.assertFalse(result.theorem_scope.applies)
        self.assertTrue(result.theorem_scope.analysis_skipped)
        self.assertEqual(result.classification, "outside_theorem_scope")
        self.assertIsNone(result.effective_polynomial)
        self.assertEqual(result.counts.total_polynomial_candidates, 0)

    def test_diagnostic_outside_theorem_scope_records_warning(self):
        result = analyze_parameter_point(
            0.2,
            -0.3,
            0.4,
            0.0,
            analyze_outside_theorem_scope=True,
        )

        self.assertFalse(result.theorem_scope.applies)
        self.assertTrue(result.theorem_scope.diagnostic_only)
        self.assertEqual(result.classification, "outside_theorem_scope")
        self.assertIn("analysis_outside_theorem_scope_is_diagnostic_only", result.warnings)

    def test_identically_zero_polynomial_classification(self):
        effective = build_effective_polynomial(np.zeros(15))
        with patch("src.workflow.effective_polynomial_from_parameters", return_value=effective):
            result = analyze_parameter_point(*GENERIC_PARAMS)

        self.assertEqual(result.classification, "identically_zero_polynomial")
        self.assertEqual(result.counts.total_polynomial_candidates, 0)

    def test_constant_nonzero_polynomial_classification(self):
        effective = build_effective_polynomial([2.0, 1e-16])
        with patch("src.workflow.effective_polynomial_from_parameters", return_value=effective):
            result = analyze_parameter_point(*GENERIC_PARAMS)

        self.assertEqual(result.classification, "constant_nonzero_polynomial")
        self.assertEqual(result.counts.total_polynomial_candidates, 0)

    def test_classification_only_domain_invalid_candidates(self):
        classification = classify_parameter_point(
            _scope(),
            build_effective_polynomial([1.0, 1.0]),
            _counts(total_polynomial_candidates=3, domain_valid_candidates=0, domain_invalid_candidates=3),
            None,
        )

        self.assertEqual(classification, "only_domain_invalid_candidates")

    def test_classification_genuine_roots_nonadmissible(self):
        classification = classify_parameter_point(
            _scope(),
            build_effective_polynomial([1.0, 1.0]),
            _counts(
                mathematically_genuine_roots=2,
                physically_admissible_roots=0,
                marginal_spatial_decay_roots=0,
                physically_nonadmissible_roots=2,
            ),
            None,
        )

        self.assertEqual(classification, "genuine_roots_nonadmissible")

    def test_classification_marginal_spatial_decay(self):
        classification = classify_parameter_point(
            _scope(),
            build_effective_polynomial([1.0, 1.0]),
            _counts(
                mathematically_genuine_roots=2,
                physically_admissible_roots=0,
                marginal_spatial_decay_roots=1,
                physically_nonadmissible_roots=1,
            ),
            None,
        )

        self.assertEqual(classification, "marginal_spatial_decay")

    def test_classification_numerically_unresolved(self):
        classification = classify_parameter_point(
            _scope(),
            build_effective_polynomial([1.0, 1.0]),
            _counts(ambiguous_mathematical_candidates=1),
            None,
        )

        self.assertEqual(classification, "numerically_unresolved")

    def test_workflow_counts_match_lower_level_outputs(self):
        result = analyze_parameter_point(*GENERIC_PARAMS)
        effective = result.effective_polynomial
        roots = candidate_roots(effective)
        verification = verify_candidate_roots(
            roots,
            *GENERIC_PARAMS,
            domain_results=result.domain_results,
        )
        physical = assess_physical_admissibility_many(verification.results)

        self.assertEqual(result.counts.total_polynomial_candidates, roots.size)
        self.assertEqual(
            result.counts.domain_valid_candidates,
            verification.summary.domain_valid_candidates,
        )
        self.assertEqual(
            result.counts.mathematically_genuine_roots,
            verification.summary.mathematically_genuine_roots,
        )
        self.assertEqual(
            result.counts.physically_admissible_roots,
            physical.summary.physically_admissible_roots,
        )


if __name__ == "__main__":
    unittest.main()
