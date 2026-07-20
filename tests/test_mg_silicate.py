import unittest

import numpy as np

from systems.mg.silicate import SilicateMgSystem, SilicateWeatheringParams


class SilicateMgSystemTests(unittest.TestCase):
    def setUp(self):
        self.system = SilicateMgSystem(basin="changjiang")

    def test_hu_example_has_physical_retained_fraction(self):
        result = self.system.calculate_weathering_flux(-0.10)

        self.assertTrue(result.success)
        self.assertEqual(result.model_status, "within_hu_calibration")
        self.assertGreater(result.f_Mg, 0.5)
        self.assertLess(result.f_Mg, 0.6)
        self.assertAlmostEqual(result.Mg_loss_fraction, 1.0 - result.f_Mg)

    def test_residue_heavier_than_protolith_means_more_mg_loss(self):
        weak = self.system.calculate_weathering_flux(-0.15)
        strong = self.system.calculate_weathering_flux(0.50)

        self.assertLess(strong.f_Mg, weak.f_Mg)
        self.assertGreater(strong.Mg_loss_fraction, weak.Mg_loss_fraction)
        self.assertEqual(strong.model_status, "extrapolated_advanced_weathering")

    def test_lighter_residue_is_not_clipped_to_unphysical_value(self):
        result = self.system.calculate_weathering_flux(-0.40)

        self.assertTrue(result.success)
        self.assertTrue(np.isnan(result.f_Mg))
        self.assertTrue(np.isnan(result.F_silicate))
        self.assertEqual(result.model_status, "incompatible_with_assumed_protolith")

    def test_released_flux_and_river_mass_balance(self):
        result = self.system.calculate_weathering_flux(-0.10)

        self.assertGreaterEqual(result.d26Mg_silicate, -0.60)
        self.assertLessEqual(result.d26Mg_silicate, -0.25)
        self.assertGreaterEqual(result.silicate_flux_fraction, 0.0)
        self.assertLessEqual(result.silicate_flux_fraction, 1.0)
        self.assertLess(result.mass_balance_check, 1e-12)
        self.assertAlmostEqual(
            result.F_silicate + result.F_carbonate,
            self.system.params.F_river_total,
        )

    def test_absolute_flux_scales_with_river_flux_prior(self):
        baseline = self.system.calculate_weathering_flux(-0.10)
        doubled = SilicateMgSystem(
            params=SilicateWeatheringParams(F_river_total=60e10)
        ).calculate_weathering_flux(-0.10)

        self.assertAlmostEqual(doubled.F_silicate, 2.0 * baseline.F_silicate)
        self.assertAlmostEqual(
            doubled.silicate_flux_fraction,
            baseline.silicate_flux_fraction,
        )

    def test_monte_carlo_is_reproducible_and_bounded(self):
        first = self.system.monte_carlo_analysis(
            -0.10,
            d26Mg_clay_std=0.02,
            n_iterations=500,
            random_seed=7,
        )
        second = self.system.monte_carlo_analysis(
            -0.10,
            d26Mg_clay_std=0.02,
            n_iterations=500,
            random_seed=7,
        )

        self.assertAlmostEqual(
            first["F_silicate"]["median"],
            second["F_silicate"]["median"],
        )
        fractions = first["silicate_flux_fraction"]["samples"]
        self.assertTrue(np.all((fractions >= 0.0) & (fractions <= 1.0)))

    def test_profile_multiplier_cancels_absolute_flux_scale(self):
        values = np.array([-0.15, -0.10, 0.00, 0.20])
        uncertainty = np.full(values.size, 0.02)
        baseline = self.system.monte_carlo_profile(
            values,
            uncertainty,
            baseline_count=2,
            n_iterations=500,
            random_seed=9,
        )
        doubled = SilicateMgSystem(
            params=SilicateWeatheringParams(F_river_total=60e10)
        ).monte_carlo_profile(
            values,
            uncertainty,
            baseline_count=2,
            n_iterations=500,
            random_seed=9,
        )

        np.testing.assert_allclose(
            baseline["flux_multiplier_median"],
            doubled["flux_multiplier_median"],
        )
        np.testing.assert_allclose(
            baseline["weathering_flux_multiplier_median"],
            doubled["weathering_flux_multiplier_median"],
        )
        np.testing.assert_allclose(
            doubled["conditional_F_silicate_median"],
            2.0 * baseline["conditional_F_silicate_median"],
        )
        self.assertTrue(baseline["absolute_flux_scale_cancels"])
        self.assertTrue(
            baseline["weathering_flux_assumes_constant_rock_supply_and_Mg"]
        )

    def test_release_profile_and_change_point_do_not_depend_on_baseline(self):
        values = np.array([-0.20, -0.18, -0.19, -0.17, 0.08, 0.10, 0.12, 0.09])
        uncertainty = np.full(values.size, 0.01)
        shallow_baseline = self.system.monte_carlo_profile(
            values,
            uncertainty,
            baseline_count=2,
            baseline_position="start",
            n_iterations=500,
            random_seed=11,
        )
        deep_baseline = self.system.monte_carlo_profile(
            values,
            uncertainty,
            baseline_count=3,
            baseline_position="end",
            n_iterations=500,
            random_seed=11,
        )

        np.testing.assert_allclose(
            shallow_baseline["Mg_release_fraction_median"],
            deep_baseline["Mg_release_fraction_median"],
        )
        self.assertEqual(
            shallow_baseline["change_point"]["after_sample_best"],
            deep_baseline["change_point"]["after_sample_best"],
        )
        self.assertEqual(
            shallow_baseline["change_point"]["after_sample_best"],
            4,
        )
        self.assertGreater(
            shallow_baseline["change_point"]["deep_to_shallow_ratio"]["median"],
            1.5,
        )

    def test_protolith_prior_is_sampled_and_can_be_fixed(self):
        uncertain = self.system.monte_carlo_analysis(
            -0.35,
            d26Mg_clay_std=0.0,
            n_iterations=1000,
            random_seed=13,
        )
        fixed = SilicateMgSystem(
            params=SilicateWeatheringParams(
                d26Mg_UCC=-0.25,
                d26Mg_UCC_range=None,
            )
        ).monte_carlo_analysis(
            -0.35,
            d26Mg_clay_std=0.0,
            n_iterations=1000,
            random_seed=13,
        )

        protolith = uncertain["d26Mg_protolith"]["samples"]
        self.assertGreaterEqual(protolith.min(), -0.45)
        self.assertLessEqual(protolith.max(), -0.25)
        self.assertGreater(uncertain["valid_weathering_fraction"], 0.4)
        self.assertEqual(fixed["valid_weathering_fraction"], 0.0)


if __name__ == "__main__":
    unittest.main()
