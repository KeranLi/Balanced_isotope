import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import image as mpimg


PROJECT_ROOT = Path(__file__).resolve().parents[1]


class MgCliProfileTests(unittest.TestCase):
    def test_siliciclastic_file_preserves_depth_order(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "profile.csv"
            output_path = Path(tmp_dir) / "result.csv"
            plot_path = Path(tmp_dir) / "profile.png"
            flux_plot_path = Path(tmp_dir) / "conditional_flux.png"
            pd.DataFrame(
                {
                    "delta_25_Mg_iso": [-0.0521, 0.0, 0.2605],
                    "delta_25_Mg_iso_2sd": [0.02, 0.02, 0.02],
                    "delta_26_Mg_iso": [-0.10, 0.0, 0.50],
                    "delta_26_Mg_iso_2sd": [0.04, 0.04, 0.04],
                }
            ).to_csv(input_path, index=False)

            completed = subprocess.run(
                [
                    sys.executable,
                    "cli.py",
                    "mg",
                    "--component-type",
                    "siliciclastic",
                    "--file",
                    str(input_path),
                    "--mg-monte-carlo",
                    "100",
                    "--output",
                    str(output_path),
                    "--plot-output",
                    str(plot_path),
                    "--flux-plot-output",
                    str(flux_plot_path),
                ],
                cwd=PROJECT_ROOT,
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = pd.read_csv(output_path)
            np.testing.assert_array_equal(result["sample_index"], [1, 2, 3])
            np.testing.assert_allclose(result["relative_depth"], [0.0, 0.5, 1.0])
            self.assertTrue(result["mass_dependent_qc"].all())
            self.assertTrue(result["absolute_flux_is_conditional"].all())
            self.assertTrue(result["protolith_is_shared_across_profile"].all())
            self.assertTrue(
                (result["delta26Mg_protolith_prior_distribution"] == "uniform").all()
            )
            np.testing.assert_allclose(
                result["delta26Mg_protolith_prior_low"], -0.45
            )
            np.testing.assert_allclose(
                result["delta26Mg_protolith_prior_high"], -0.25
            )
            self.assertTrue(result["F_silicate_mc_median_mol_yr"].notna().all())
            self.assertTrue(
                result["Mg_release_fraction_mc_median"].between(0.0, 1.0).all()
            )
            self.assertTrue(
                result[
                    "conditional_F_silicate_profile_mc_median_mol_yr"
                ].notna().all()
            )
            self.assertIn("automatic_change_point_after_sample", result.columns)
            self.assertIn(
                "deep_to_shallow_weathering_flux_ratio_median",
                result.columns,
            )
            self.assertNotIn(
                "relative_weathering_flux_multiplier_median",
                result.columns,
            )
            image = mpimg.imread(plot_path)
            self.assertGreater(image.shape[0], image.shape[1])
            flux_image = mpimg.imread(flux_plot_path)
            self.assertGreater(flux_image.shape[0], flux_image.shape[1])


if __name__ == "__main__":
    unittest.main()
