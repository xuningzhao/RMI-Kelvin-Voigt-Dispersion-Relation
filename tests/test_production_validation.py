from __future__ import annotations

import unittest

from src.production_validation import conjugate_summary, spectra_match, verify_residuals
from src.workflow import analyze_parameter_point


class ProductionValidationTests(unittest.TestCase):
    def test_conjugate_pairing_reports_real_pairs_and_unmatched(self):
        summary = conjugate_summary([-1 + 0j, -2 + 3j, -2 - 3j, 1 + 2j])
        self.assertEqual((summary.real_count, summary.pair_count), (1, 1))
        self.assertEqual(summary.unmatched, (1 + 2j,))

    def test_residual_verification_accepts_production_roots(self):
        result = analyze_parameter_point(0.2, -0.3, 0.4, 1.7)
        roots = result.verification_run.genuine_roots
        residuals = verify_residuals(roots, 0.2, -0.3, 0.4, 1.7)
        self.assertEqual(len(residuals), len(roots))

    def test_material_exchange_spectra_match(self):
        first = analyze_parameter_point(0.2, -0.3, 0.4, 1.7)
        second = analyze_parameter_point(-0.2, 0.3, -0.4, 1.7)
        self.assertTrue(spectra_match(first.verification_run.genuine_roots, second.verification_run.genuine_roots))


if __name__ == "__main__":
    unittest.main()
