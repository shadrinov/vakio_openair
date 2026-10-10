import unittest

import numpy as np

from scripts.axisymmetric_vakio import solve


class AxisymmetricTests(unittest.TestCase):
    def test_equilibrium_with_equal_boundary_temperatures(self):
        result = solve(0.036, step_m=0.01, indoor_c=25, outdoor_c=25)
        np.testing.assert_allclose(result.temperatures_c, 25, atol=1e-7)
        np.testing.assert_allclose(result.room_surface_c, 25, atol=1e-7)

    def test_energy_balance_and_physical_temperature_bounds(self):
        result = solve(0.036, step_m=0.01)
        self.assertLess(result.relative_heat_imbalance, 1e-6)
        self.assertGreaterEqual(result.temperatures_c.min(), -20)
        self.assertLessEqual(result.temperatures_c.max(), 25)

    def test_surface_result_converges_under_mesh_refinement(self):
        coarse = solve(0.036, step_m=0.005)
        fine = solve(0.036, step_m=0.0025)
        self.assertLess(abs(coarse.temperature_at_radius(0.0895) - fine.temperature_at_radius(0.0895)), 0.05)
