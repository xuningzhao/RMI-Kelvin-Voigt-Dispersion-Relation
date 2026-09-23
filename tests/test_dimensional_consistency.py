"""Independent dimensional-to-nondimensional consistency tests.

Multiplication of the dimensional equation by ``rhoT / k**2`` gives the
pre-contrast nondimensional equation with constant 4.  Substituting the
contrast fractions contributes a factor 2 to both reciprocal terms, after
which the adopted final equation is divided by 2 and has constant 2.
Therefore the exact proportionality tested here is

    D_star_Lambda = (rhoT / (2*k**2)) * original_dimensional_D.

This factor is derived from the equation normalization, not fitted from data.
"""

from __future__ import annotations

import unittest

import numpy as np

from src.dispersion import (
    D_star_Ek,
    D_star_Lambda,
    nondimensional_groups,
    original_dimensional_D,
)


class DimensionalConsistencyTests(unittest.TestCase):
    def test_conversion_helper_matches_definitions(self):
        dimensional = (0.4 + 0.7j, 1.3, 2.0, 5.0, 0.8, 1.7, 3.0, 0.25)
        s, arho, amu, ag, lam, ek = nondimensional_groups(*dimensional)
        gamma, k, rho1, rho2, mu1, mu2, g1, g2 = dimensional
        rho_total = rho1 + rho2
        mu_total = mu1 + mu2
        g_total = g1 + g2

        self.assertAlmostEqual(s, rho_total * gamma / (mu_total * k**2))
        self.assertAlmostEqual(arho, (rho1 - rho2) / rho_total)
        self.assertAlmostEqual(amu, (mu1 - mu2) / mu_total)
        self.assertAlmostEqual(ag, (g1 - g2) / g_total)
        self.assertAlmostEqual(lam, rho_total * g_total / (mu_total**2 * k**2))
        self.assertAlmostEqual(lam, ((1.0 - ek) / (1.0 + ek)) ** 2)

    def test_1000_random_dimensional_cases(self):
        rng = np.random.default_rng(8675309)
        count = 1000

        # Log-uniform sampling supplies several orders of magnitude of density,
        # viscosity, modulus, and wavenumber ratios.  The first 100 cases in
        # each modulus array force one medium close to the zero-modulus limit.
        rho1 = 10.0 ** rng.uniform(-2.0, 2.0, count)
        rho2 = 10.0 ** rng.uniform(-2.0, 2.0, count)
        mu1 = 10.0 ** rng.uniform(-2.0, 2.0, count)
        mu2 = 10.0 ** rng.uniform(-2.0, 2.0, count)
        g1 = 10.0 ** rng.uniform(-5.0, 2.0, count)
        g2 = 10.0 ** rng.uniform(-5.0, 2.0, count)
        g1[:100] = 10.0 ** rng.uniform(-14.0, -10.0, 100)
        g2[100:200] = 10.0 ** rng.uniform(-14.0, -10.0, 100)
        k = 10.0 ** rng.uniform(-2.0, 2.0, count)

        gamma_scale = 10.0 ** rng.uniform(-2.0, 2.0, count)
        gamma = gamma_scale * (
            rng.uniform(-2.0, 2.0, count)
            + 1j * rng.uniform(0.05, 2.0, count)
        )

        dimensional = original_dimensional_D(
            gamma, k, rho1, rho2, mu1, mu2, g1, g2
        )
        s, arho, amu, ag, lam, ek = nondimensional_groups(
            gamma, k, rho1, rho2, mu1, mu2, g1, g2
        )
        nondimensional = D_star_Lambda(s, arho, amu, ag, lam)

        # Exact normalization derived in the module docstring above.
        scaling = (rho1 + rho2) / (2.0 * k**2)
        np.testing.assert_allclose(
            nondimensional,
            scaling * dimensional,
            rtol=2e-11,
            atol=2e-11,
        )

        # The bounded coordinate must describe the same evaluations.  Near
        # Ek = -1, converting Lambda -> Ek -> Lambda is ill-conditioned, so a
        # small floating-point round-trip tolerance is necessary.
        np.testing.assert_allclose(
            D_star_Ek(s, arho, amu, ag, ek),
            nondimensional,
            rtol=1e-12,
            atol=1e-12,
        )

        # Confirm that the sample genuinely exercises the requested regimes.
        self.assertLess(np.min(g1), 1e-10)
        self.assertLess(np.min(g2), 1e-10)
        self.assertGreater(np.max(rho1 / rho2), 1e3)
        self.assertGreater(np.max(rho2 / rho1), 1e3)
        self.assertGreater(np.max(mu1 / mu2), 1e3)
        self.assertGreater(np.max(mu2 / mu1), 1e3)
        self.assertGreater(np.max(k) / np.min(k), 1e3)
        self.assertTrue(np.all(np.imag(gamma) > 0.0))

    def test_scalar_dimensional_scaling_and_zero_individual_modulus(self):
        dimensional_inputs = (0.3 + 0.9j, 0.7, 4.0, 0.25, 1.5, 0.08, 0.0, 2.2)
        dimensional = original_dimensional_D(*dimensional_inputs)
        groups = nondimensional_groups(*dimensional_inputs)
        rho_total = dimensional_inputs[2] + dimensional_inputs[3]
        k = dimensional_inputs[1]

        expected = rho_total * dimensional / (2.0 * k**2)
        actual = D_star_Lambda(*groups[:5])
        self.assertAlmostEqual(actual.real, expected.real, places=12)
        self.assertAlmostEqual(actual.imag, expected.imag, places=12)

    def test_conversion_rejects_undefined_elastic_contrast(self):
        with self.assertRaisesRegex(ValueError, "AG is undefined"):
            nondimensional_groups(1.0j, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0)


if __name__ == "__main__":
    unittest.main()
