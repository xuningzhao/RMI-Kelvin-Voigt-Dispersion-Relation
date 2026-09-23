"""Tests for the nondimensional dispersion-function evaluators."""

from __future__ import annotations

import unittest
import warnings

import numpy as np

from src.dispersion import D_star_Ek, D_star_Lambda


def _direct_reference(s, arho, amu, ag, lam):
    """Independent transcription of the compact formula."""

    eta1 = (1.0 + amu) + (lam / s) * (1.0 + ag)
    eta2 = (1.0 - amu) + (lam / s) * (1.0 - ag)
    q1 = np.sqrt(1.0 + ((1.0 + arho) * s) / eta1)
    q2 = np.sqrt(1.0 + ((1.0 - arho) * s) / eta2)
    return s * (1.0 / (eta1 + eta2 * q2) + 1.0 / (eta2 + eta1 * q1)) + 2.0


def _cleared_reference(s, arho, amu, ag, lam):
    """Independent equation (5), used to check normalization and factors."""

    b1 = s * (1.0 + amu) + lam * (1.0 + ag)
    b2 = s * (1.0 - amu) + lam * (1.0 - ag)
    q1 = np.sqrt(1.0 + ((1.0 + arho) * s**2) / b1)
    q2 = np.sqrt(1.0 + ((1.0 - arho) * s**2) / b2)
    return s**2 * (1.0 / (b1 + b2 * q2) + 1.0 / (b2 + b1 * q1)) + 2.0


class DispersionTests(unittest.TestCase):
    def test_scalar_complex_input(self):
        s = -0.35 + 0.8j
        parameters = (0.2, -0.3, 0.45, 1.7)
        actual = D_star_Lambda(s, *parameters)
        expected = _direct_reference(s, *parameters)

        self.assertIsInstance(actual, complex)
        self.assertAlmostEqual(actual.real, expected.real, places=14)
        self.assertAlmostEqual(actual.imag, expected.imag, places=14)

    def test_numpy_array_input_matches_scalar_evaluation(self):
        s_values = np.array([-1.2 + 0.1j, -0.4 + 0.7j, 0.3 + 1.1j])
        array_result = D_star_Lambda(s_values, 0.1, -0.2, 0.3, 0.9)
        scalar_result = np.array(
            [D_star_Lambda(value, 0.1, -0.2, 0.3, 0.9) for value in s_values]
        )

        self.assertIsInstance(array_result, np.ndarray)
        self.assertEqual(array_result.shape, s_values.shape)
        np.testing.assert_allclose(array_result, scalar_result, rtol=1e-14, atol=1e-14)

    def test_random_parameters_match_cleared_equation(self):
        rng = np.random.default_rng(20260706)
        for _ in range(100):
            s = rng.uniform(-2.0, 2.0) + 1j * rng.uniform(0.1, 2.0)
            arho, amu, ag = rng.uniform(-0.9, 0.9, size=3)
            lam = rng.uniform(0.0, 8.0)

            actual = D_star_Lambda(s, arho, amu, ag, lam)
            expected = _cleared_reference(s, arho, amu, ag, lam)
            np.testing.assert_allclose(actual, expected, rtol=2e-13, atol=2e-13)

    def test_Ek_to_Lambda_consistency(self):
        rng = np.random.default_rng(1729)
        s = rng.uniform(-2.0, 2.0, size=128) + 1j * rng.uniform(0.1, 2.0, size=128)
        arho = rng.uniform(-0.95, 0.95, size=128)
        amu = rng.uniform(-0.95, 0.95, size=128)
        ag = rng.uniform(-0.95, 0.95, size=128)
        ek = rng.uniform(-0.95, 1.0, size=128)
        lam = ((1.0 - ek) / (1.0 + ek)) ** 2

        from_ek = D_star_Ek(s, arho, amu, ag, ek)
        from_lambda = D_star_Lambda(s, arho, amu, ag, lam)
        np.testing.assert_allclose(from_ek, from_lambda, rtol=0.0, atol=0.0)

    def test_Ek_endpoint_at_one_matches_zero_Lambda(self):
        s = np.array([-0.5 + 0.2j, 0.7 + 0.9j])
        np.testing.assert_allclose(
            D_star_Ek(s, 0.2, -0.1, 0.4, 1.0),
            D_star_Lambda(s, 0.2, -0.1, 0.4, 0.0),
            rtol=0.0,
            atol=0.0,
        )

    def test_medium_exchange_symmetry(self):
        s = np.array([-0.8 + 0.4j, 0.2 + 1.3j])
        original = D_star_Lambda(s, 0.6, -0.25, 0.15, 2.3)
        exchanged = D_star_Lambda(s, -0.6, 0.25, -0.15, 2.3)
        np.testing.assert_allclose(original, exchanged, rtol=1e-14, atol=1e-14)

    def test_zero_s_is_rejected_for_scalar_and_array(self):
        with self.assertRaisesRegex(ValueError, "s = 0 is singular"):
            D_star_Lambda(0.0, 0.0, 0.0, 0.0, 1.0)

        with self.assertRaisesRegex(ValueError, "s = 0 is singular"):
            D_star_Lambda(np.array([1.0j, 0.0]), 0.0, 0.0, 0.0, 1.0)

    def test_small_nonzero_s_is_not_silently_divided_by_zero(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            result = D_star_Lambda(1e-14 + 2e-14j, 0.2, -0.1, 0.3, 0.8)

        self.assertTrue(np.isfinite(result.real))
        self.assertTrue(np.isfinite(result.imag))

    def test_invalid_parameter_domains_are_rejected(self):
        with self.assertRaises(ValueError):
            D_star_Lambda(1.0j, 1.01, 0.0, 0.0, 1.0)
        with self.assertRaises(ValueError):
            D_star_Lambda(1.0j, 0.0, 0.0, 0.0, -1.0)
        with self.assertRaises(ValueError):
            D_star_Ek(1.0j, 0.0, 0.0, 0.0, -1.0)


if __name__ == "__main__":
    unittest.main()

