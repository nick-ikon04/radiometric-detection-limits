import unittest

from radiometric_detection_limits import calc_peak_roi, calc_radiometry, normal_quantile


class RadiometryTests(unittest.TestCase):
    def test_standard_normal_quantiles(self):
        self.assertAlmostEqual(normal_quantile(0.95), 1.6448536269514722)
        self.assertAlmostEqual(normal_quantile(0.975), 1.9599639845400536)

    def test_surface_smear_example(self):
        result = calc_radiometry(
            background_counts=1652,
            background_time_s=1800,
            gross_counts=2861,
            gross_time_s=60,
            efficiency=0.47,
            regulatory_level=200,
        )
        self.assertAlmostEqual(result.activity_bq, 99.50118203309692)
        self.assertTrue(result.signal_detected)
        self.assertTrue(result.method_suitable)

    def test_radioactive_waste_example(self):
        result = calc_radiometry(
            background_counts=256,
            background_time_s=7200,
            gross_counts=779,
            gross_time_s=3600,
            efficiency=0.018,
            mass_kg=0.750,
            beta=0.001,
            confidence=0.998,
            regulatory_level=50,
        )
        self.assertAlmostEqual(result.activity_bq, 10.046296296296296)
        self.assertAlmostEqual(result.specific_activity_bq_per_kg, 13.39506172839506)
        self.assertTrue(result.signal_detected)
        self.assertTrue(result.method_suitable)

    def test_be7_roi_example(self):
        counts = dict(
            zip(
                range(378, 398),
                [22, 19, 21, 12, 15, 11, 13, 17, 20, 25, 26, 36, 20, 20, 14, 11, 13, 18, 18, 19],
            )
        )
        result = calc_peak_roi(
            counts,
            (384, 390),
            ((380, 383), (391, 394)),
            50_000,
            0.01315,
        )
        self.assertAlmostEqual(result.gross_peak_counts, 157.0)
        self.assertAlmostEqual(result.background_under_peak_counts, 102.375)
        self.assertAlmostEqual(result.net_peak_counts, 54.625)
        self.assertAlmostEqual(result.activity_bq, 0.08307984790874524)
        self.assertTrue(result.signal_detected)

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(ValueError):
            calc_radiometry(1, 0, 1, 10, 0.5)
        with self.assertRaises(ValueError):
            calc_peak_roi({1: 2}, (1, 2), ((2, 3),), 10, 0.5)


if __name__ == "__main__":
    unittest.main()

