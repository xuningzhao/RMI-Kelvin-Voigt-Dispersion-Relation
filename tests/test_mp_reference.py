"""Reference-solver invariants and manuscript air/hydrogel regression."""
import unittest

import mpmath as mp

from src.mp_reference import solve_reference_point
from src.dispersion import D_star_Lambda


class ReferenceSolverTests(unittest.TestCase):
    def test_air_hydrogel_oscillatory_reference_and_precision_convergence(self):
        inputs = ("-0.9974", "-0.999998", "-1", "0")
        low = solve_reference_point(*inputs, decimal_precision=60)
        high = solve_reference_point(*inputs, decimal_precision=80)
        self.assertEqual(len(low.accepted), 2)
        self.assertEqual(len(high.accepted), 2)
        for s, qp, qm, residual in high.accepted:
            self.assertLess(abs(s.real - (-0.45678649147107298)), 1e-12)
            self.assertLess(abs(abs(s.imag) - 0.83959400219637381), 1e-12)
            self.assertLess(abs(D_star_Lambda(s, -0.9974, -0.999998, -1, 1)), 1e-8)
            self.assertGreater(qp.real, 0)
            self.assertGreater(qm.real, 0)
            self.assertLess(residual, 1e-40)
            self.assertLess(min(abs(s - row[0]) for row in low.accepted), 1e-12)
            self.assertLess(min(abs(s.conjugate() - row[0]) for row in high.accepted), 1e-12)

    def test_medium_exchange_preserves_accepted_spectrum(self):
        original = solve_reference_point("0.2", "-0.3", "0.4", "0.1")
        exchanged = solve_reference_point("-0.2", "0.3", "-0.4", "0.1")
        self.assertGreater(len(original.accepted), 0)
        self.assertEqual(len(original.accepted), len(exchanged.accepted))
        for s, qp, qm, _ in original.accepted:
            other = min(exchanged.accepted, key=lambda row: abs(row[0] - s))
            self.assertLess(abs(other[0] - s), 1e-12)
            self.assertLess(abs(other[1] - qm), 1e-12)
            self.assertLess(abs(other[2] - qp), 1e-12)

    def test_precision_context_restored(self):
        with mp.workdps(35):
            solve_reference_point("0", "0", "-1", "0.9")
            self.assertEqual(mp.mp.dps, 35)

    def test_outside_scope_and_invalid_precision_rejected(self):
        for values in (("1", "0", "0", "0"), ("0", "0", "0", "1"),
                       ("0", "0", "0", "-1"), ("nan", "0", "0", "0")):
            with self.subTest(values=values), self.assertRaises(ValueError):
                solve_reference_point(*values)
        for precision in (0, 20, 80.5, True):
            with self.subTest(precision=precision), self.assertRaises(ValueError):
                solve_reference_point("0", "0", "0", "0", decimal_precision=precision)
