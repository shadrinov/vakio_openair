import unittest

from dewpoint.model import BLOCK_THICKNESS_MM, calculate_wall, dew_point


class WallModelTests(unittest.TestCase):
    def test_saturated_air_has_dew_point_equal_to_air_temperature(self):
        self.assertAlmostEqual(dew_point(20, 100), 20)

    def test_winter_example_crosses_inside_block(self):
        result = calculate_wall(20, -10, 50)
        self.assertAlmostEqual(result.dew_point_c, 9.27, places=2)
        self.assertIsNotNone(result.isotherm_mm)
        self.assertLess(result.isotherm_mm, BLOCK_THICKNESS_MM)
        self.assertGreater(result.temperatures_c[0], result.dew_point_c)

    def test_warm_outside_has_no_crossing(self):
        self.assertIsNone(calculate_wall(20, 10, 50).isotherm_mm)

    def test_invalid_humidity_is_rejected(self):
        with self.assertRaises(ValueError):
            dew_point(20, 0)


if __name__ == "__main__":
    unittest.main()
