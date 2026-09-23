"""Tests for the P14 polynomial-candidate infrastructure."""

from __future__ import annotations

from pathlib import Path
import unittest

import numpy as np
import sympy as sp
from sympy.parsing.mathematica import parse_mathematica

from src.polynomial import (
    build_effective_polynomial,
    candidate_roots,
    effective_polynomial_from_parameters,
    evaluate_coefficients,
    merge_duplicate_roots,
    polynomial_value,
)
from src.parameter_sweep import lambda_from_ek


ROOT = Path(__file__).resolve().parents[1]
COEFF_EXPORT = ROOT / "symbolic" / "exports" / "canonical" / "P14_coefficients.wl"
POLY_EXPORT = ROOT / "symbolic" / "exports" / "canonical" / "P14_polynomial.wl"


REFERENCE_CASES = (
    (0.2, -0.3, 0.4, 1.7, 14, "generic-degree-14"),
    (0.3, -0.3, 0.2, 0.8, 13, "reduced-leading-degree"),
    (-0.65, 0.15, -0.5, 2.25, 14, "nonsymmetric-degree-14"),
)

NEGATIVE_EK_REGRESSION_CASES = (
    (-0.8, -0.8, -0.44209312324515126 + 8.451250773470967j),
    (-0.64, -0.64, -0.446661767207 + 4.13308125598j),
    (0.64, -0.8, -2.74832226154 - 11.3650958193j),
    (0.64, -0.64, -1.62000870959 + 5.30724165312j),
)


def _mathematica_coefficients():
    text = COEFF_EXPORT.read_text().replace("Lambda", "Lam")
    return parse_mathematica(text)


def _mathematica_polynomial():
    text = POLY_EXPORT.read_text().replace("Lambda", "Lam")
    return parse_mathematica(text)


def _reference_coefficients(expressions, Arho, Amu, AG, Lambda):
    symbols = sp.symbols("Arho Amu AG Lam")
    substitutions = dict(zip(symbols, (Arho, Amu, AG, Lambda), strict=True))
    return np.asarray(
        [complex(sp.N(expr.subs(substitutions), 40)) for expr in expressions],
        dtype=np.complex128,
    )


class PolynomialInfrastructureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mathematica_coefficients = _mathematica_coefficients()
        cls.mathematica_polynomial = _mathematica_polynomial()

    def test_coefficients_match_frozen_mathematica_export(self):
        for Arho, Amu, AG, Lambda, expected_degree, name in REFERENCE_CASES:
            with self.subTest(name=name):
                python_coeffs = evaluate_coefficients(Arho, Amu, AG, Lambda)
                reference_coeffs = _reference_coefficients(
                    self.mathematica_coefficients, Arho, Amu, AG, Lambda
                )

                self.assertEqual(python_coeffs.shape, (15,))
                np.testing.assert_allclose(
                    python_coeffs, reference_coeffs, rtol=2e-13, atol=2e-10
                )

                effective = build_effective_polynomial(
                    python_coeffs, relative_tolerance=1e-12
                )
                self.assertEqual(effective.effective_degree, expected_degree)

    def test_polynomial_value_matches_frozen_mathematica_polynomial(self):
        s_value = 0.37 - 0.91j
        Arho_symbol, Amu_symbol, AG_symbol, Lam_symbol, s_symbol = sp.symbols(
            "Arho Amu AG Lam s"
        )

        for Arho, Amu, AG, Lambda, _expected_degree, name in REFERENCE_CASES:
            with self.subTest(name=name):
                python_value = polynomial_value(
                    s_value, evaluate_coefficients(Arho, Amu, AG, Lambda)
                )
                reference = complex(
                    sp.N(
                        self.mathematica_polynomial.subs(
                            {
                                Arho_symbol: Arho,
                                Amu_symbol: Amu,
                                AG_symbol: AG,
                                Lam_symbol: Lambda,
                                s_symbol: s_value,
                            }
                        ),
                        40,
                    )
                )

                self.assertAlmostEqual(python_value.real, reference.real, places=9)
                self.assertAlmostEqual(python_value.imag, reference.imag, places=9)

    def test_effective_polynomial_classifies_constant_and_zero_cases(self):
        zero = build_effective_polynomial(np.zeros(15), relative_tolerance=1e-12)
        self.assertEqual(zero.status, "identically_zero")
        self.assertIsNone(zero.effective_degree)

        constant = build_effective_polynomial([3.0, 1e-16], relative_tolerance=1e-12)
        self.assertEqual(constant.status, "constant_nonzero")
        self.assertEqual(constant.effective_degree, 0)
        self.assertEqual(candidate_roots(constant).size, 0)

    def test_candidate_roots_have_small_polynomial_residuals(self):
        effective = effective_polynomial_from_parameters(0.2, -0.3, 0.4, 1.7)
        roots = candidate_roots(effective)

        self.assertEqual(effective.status, "generic_degree")
        self.assertEqual(roots.size, 14)

        residuals = [
            abs(polynomial_value(root, effective.coefficients_ascending))
            for root in roots
        ]
        self.assertLess(max(residuals) / effective.coefficient_scale, 1e-5)

    def test_high_lambda_negative_Ek_preserves_full_structural_degree(self):
        for AG, Ek, reference_root in NEGATIVE_EK_REGRESSION_CASES:
            with self.subTest(AG=AG, Ek=Ek):
                Lambda = lambda_from_ek(Ek)
                effective = effective_polynomial_from_parameters(-0.95, -0.8, AG, Lambda)
                roots = candidate_roots(effective)
                nearest = min(roots, key=lambda root: abs(root - reference_root))

                self.assertEqual(effective.status, "generic_degree")
                self.assertEqual(effective.effective_degree, 14)
                self.assertEqual(roots.size, 14)
                self.assertLess(abs(nearest - reference_root), 1e-8)
                self.assertLess(
                    abs(polynomial_value(reference_root, effective.coefficients_ascending))
                    / effective.coefficient_scale,
                    1e-8,
                )

    def test_duplicate_root_merging_preserves_multiplicity_information(self):
        roots = np.asarray([1 + 2j, 1 + 2j + 1e-10, -0.5 + 0.25j])
        clusters = merge_duplicate_roots(roots, absolute_tolerance=1e-8)

        multiplicities = sorted(cluster.multiplicity for cluster in clusters)
        self.assertEqual(multiplicities, [1, 2])
        self.assertEqual(sum(cluster.multiplicity for cluster in clusters), roots.size)


if __name__ == "__main__":
    unittest.main()
