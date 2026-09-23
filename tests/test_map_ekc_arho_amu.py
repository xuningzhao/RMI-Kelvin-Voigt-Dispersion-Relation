import numpy as np
import unittest

from scripts.map_ekc_arho_amu import _down_crossings


class DownCrossingTests(unittest.TestCase):
    def test_down_crossings_returns_each_complex_to_noncomplex_interval(self):
        states = np.asarray([False, True, True, False, True, False, False])
        np.testing.assert_array_equal(_down_crossings(states), np.asarray([2, 4]))


    def test_down_crossings_handles_no_transition(self):
        np.testing.assert_array_equal(_down_crossings(np.ones(5, dtype=bool)), np.asarray([]))
