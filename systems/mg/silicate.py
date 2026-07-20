"""Siliciclastic Mg-isotope weathering model.

The implementation follows the two distinct calculations in Hu et al.
(2023):

1. A Rayleigh model relates the isotope composition of the weathering
   residue to the fraction of Mg retained in that residue.
2. A riverine two-endmember mass balance partitions dissolved Mg between
   silicate and carbonate weathering sources.

Absolute fluxes are conditional on ``F_river_total``. A clay isotope value
alone constrains a weathering state, not an absolute flux in mol/yr.
"""

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Tuple

import numpy as np


@dataclass
class SilicateWeatheringParams:
    """Parameters for the Hu et al. (2023) siliciclastic model."""

    # Isotope endmembers (per mil vs DSM3)
    d26Mg_UCC: float = -0.25
    d26Mg_UCC_range: Optional[Tuple[float, float]] = (-0.45, -0.25)
    d26Mg_carbonate: float = -2.0
    d26Mg_river_water: float = -1.14

    # Dissolved Mg relative to fresh silicate. Hu et al. constrain this to
    # -0.35 to -0.15 per mil for incipient-to-intermediate weathering.
    Delta_fluid_protolith: float = -0.25
    Delta_fluid_protolith_range: Tuple[float, float] = (-0.35, -0.15)

    # Total dissolved river Mg flux. This is the Changjiang prior used by
    # Hu et al.; ancient-section fluxes scale linearly with this assumption.
    F_river_total: float = 30e10

    # One-sigma uncertainties unless explicitly stated otherwise.
    uncertainty_clay: float = 0.07
    uncertainty_river: float = 0.15
    uncertainty_carbonate: float = 0.40  # Hu et al.: -2.0 +/- 0.8 (2SD)
    uncertainty_river_flux_fraction: float = 0.0

    # Hu et al. calibrate the Rayleigh model for less than 50% Mg depletion.
    calibration_min_f_Mg: float = 0.5
    weathering_stages: Dict[str, Tuple[float, float]] = field(
        default_factory=lambda: {
            "incipient": (0.8, 1.0),
            "intermediate": (0.5, 0.8),
            "advanced": (0.2, 0.5),
            "extreme": (0.0, 0.2),
        }
    )


@dataclass
class SilicateModelResult:
    """Result from one siliciclastic Mg-isotope calculation."""

    success: bool = True
    message: str = ""
    data: Dict[str, Any] = field(default_factory=dict)

    f_Mg: Optional[float] = None
    Mg_loss_fraction: Optional[float] = None
    weathering_stage: Optional[str] = None
    model_status: Optional[str] = None
    within_hu_calibration: bool = False

    d26Mg_silicate: Optional[float] = None
    d26Mg_carbonate: float = -2.0

    silicate_flux_fraction: Optional[float] = None
    F_silicate: Optional[float] = None
    F_carbonate: Optional[float] = None
    SWI: Optional[float] = None
    mass_balance_check: Optional[float] = None

    def get(self, key: str, default=None):
        return self.data.get(key, default)


class SilicateMgSystem:
    """Infer silicate weathering state and conditional dissolved Mg flux."""

    COMPONENT_TYPE = "siliciclastic"
    ELEMENT = "mg"
    NAME = "Magnesium (Siliciclastic)"

    def __init__(
        self,
        params: Optional[SilicateWeatheringParams] = None,
        basin: str = "changjiang",
    ):
        self.params = params or self._get_default_params(basin)
        self.basin = basin
        self._validate_parameters()
        self._last_result: Optional[SilicateModelResult] = None

    def _get_default_params(self, basin: str) -> SilicateWeatheringParams:
        if basin == "changjiang":
            return SilicateWeatheringParams()
        if basin == "global":
            return SilicateWeatheringParams(
                d26Mg_river_water=-1.0,
                F_river_total=100e10,
            )
        return SilicateWeatheringParams()

    def _validate_parameters(self) -> None:
        p = self.params
        if p.Delta_fluid_protolith >= 0:
            raise ValueError("Delta_fluid_protolith must be negative")
        if p.d26Mg_UCC <= -1000:
            raise ValueError("d26Mg_UCC must be greater than -1000")
        if p.F_river_total <= 0:
            raise ValueError("F_river_total must be positive")
        low, high = p.Delta_fluid_protolith_range
        if low >= high or high >= 0:
            raise ValueError("Delta_fluid_protolith_range must be increasing and negative")
        if p.d26Mg_UCC_range is not None:
            low, high = p.d26Mg_UCC_range
            if low >= high or low <= -1000:
                raise ValueError(
                    "d26Mg_UCC_range must be increasing and greater than -1000"
                )

    def _draw_protolith(
        self,
        rng: np.random.Generator,
        n_iterations: int,
    ) -> np.ndarray:
        """Draw one protolith value per realization.

        A profile uses the same one-dimensional draw array for every sample,
        so protolith uncertainty cannot create artificial stratigraphic
        variation within a realization.
        """
        if self.params.d26Mg_UCC_range is None:
            return np.full(n_iterations, self.params.d26Mg_UCC)
        return rng.uniform(*self.params.d26Mg_UCC_range, n_iterations)

    @staticmethod
    def _delta_to_relative_ratio(delta: np.ndarray) -> np.ndarray:
        return 1.0 + np.asarray(delta, dtype=float) / 1000.0

    def _raw_f_mg(
        self,
        d26Mg_clay: float,
        d26Mg_protolith: Optional[float] = None,
        delta_fluid_protolith: Optional[float] = None,
    ) -> float:
        """Invert the exact Rayleigh residual equation for retained Mg."""
        d26Mg_protolith = (
            self.params.d26Mg_UCC
            if d26Mg_protolith is None
            else d26Mg_protolith
        )
        delta_fluid_protolith = (
            self.params.Delta_fluid_protolith
            if delta_fluid_protolith is None
            else delta_fluid_protolith
        )
        alpha = 1.0 + delta_fluid_protolith / 1000.0
        if alpha <= 0 or np.isclose(alpha, 1.0):
            raise ValueError("Invalid fluid-protolith fractionation factor")

        residue_ratio = self._delta_to_relative_ratio(d26Mg_clay)
        initial_ratio = self._delta_to_relative_ratio(d26Mg_protolith)
        if residue_ratio <= 0 or initial_ratio <= 0:
            raise ValueError("Delta values must be greater than -1000 per mil")

        return float((residue_ratio / initial_ratio) ** (1.0 / (alpha - 1.0)))

    def calculate_weathering_degree(self, d26Mg_clay: float) -> Tuple[float, str]:
        """Return retained Mg fraction and weathering-stage classification.

        For a negative fluid-protolith fractionation, a residue lighter than
        the assumed protolith implies ``f_Mg > 1`` and has no physical solution
        under this model. Such values are returned as NaN rather than clipped.
        """
        f_raw = self._raw_f_mg(d26Mg_clay)
        if f_raw > 1.0 + 1e-9:
            return float("nan"), "incompatible_protolith"

        f_mg = float(np.clip(f_raw, 0.0, 1.0))
        return f_mg, self._classify_weathering_stage(f_mg)

    def _classify_weathering_stage(self, f_mg: float) -> str:
        if not np.isfinite(f_mg):
            return "incompatible_protolith"
        for stage_name, (f_min, f_max) in self.params.weathering_stages.items():
            if f_min <= f_mg <= f_max:
                return stage_name
        return "extreme"

    def calculate_released_flux_isotope(
        self,
        d26Mg_clay: float,
        f_mg: Optional[float] = None,
        d26Mg_protolith: Optional[float] = None,
        delta_fluid_protolith: Optional[float] = None,
    ) -> float:
        """Calculate isotope composition of cumulative Mg released.

        For Rayleigh removal, the accumulated product ratio is

        ``R_release = R0 * (1 - f**alpha) / (1 - f)``.

        It approaches ``alpha * R0`` at ``f = 1`` and the protolith ratio as
        Mg depletion approaches completion.
        """
        d26Mg_protolith = (
            self.params.d26Mg_UCC
            if d26Mg_protolith is None
            else d26Mg_protolith
        )
        delta_fluid_protolith = (
            self.params.Delta_fluid_protolith
            if delta_fluid_protolith is None
            else delta_fluid_protolith
        )
        if f_mg is None:
            f_mg, _ = self.calculate_weathering_degree(d26Mg_clay)
        if not np.isfinite(f_mg) or not 0 <= f_mg <= 1:
            return float("nan")

        alpha = 1.0 + delta_fluid_protolith / 1000.0
        initial_ratio = float(self._delta_to_relative_ratio(d26Mg_protolith))
        if np.isclose(f_mg, 1.0):
            released_ratio = alpha * initial_ratio
        else:
            released_ratio = initial_ratio * (1.0 - f_mg**alpha) / (1.0 - f_mg)
        return float((released_ratio - 1.0) * 1000.0)

    # Backward-compatible name. The value now represents the cumulative
    # dissolved silicate-weathering flux, not clay plus an arbitrary offset.
    def calculate_silicate_endmember(self, d26Mg_clay: float) -> float:
        return self.calculate_released_flux_isotope(d26Mg_clay)

    @staticmethod
    def _solve_flux_fraction(
        d26Mg_river: float,
        d26Mg_silicate: float,
        d26Mg_carbonate: float,
    ) -> float:
        denominator = d26Mg_silicate - d26Mg_carbonate
        if abs(denominator) < 1e-12:
            return float("nan")
        return float((d26Mg_river - d26Mg_carbonate) / denominator)

    def calculate_weathering_flux(self, d26Mg_clay: float) -> SilicateModelResult:
        """Calculate weathering state and conditional absolute Mg flux."""
        p = self.params
        f_mg, stage = self.calculate_weathering_degree(d26Mg_clay)
        f_raw = self._raw_f_mg(d26Mg_clay)

        if not np.isfinite(f_mg):
            status = "incompatible_with_assumed_protolith"
            message = (
                "Observed residue is lighter than the assumed protolith; "
                "the negative-fractionation Rayleigh model has no physical solution."
            )
        elif f_mg < p.calibration_min_f_Mg:
            status = "extrapolated_advanced_weathering"
            message = (
                "Rayleigh result is outside the Hu et al. (2023) "
                "incipient-to-intermediate calibration range."
            )
        else:
            status = "within_hu_calibration"
            message = ""

        d26Mg_silicate = self.calculate_released_flux_isotope(
            d26Mg_clay,
            f_mg=f_mg,
        )
        flux_fraction = self._solve_flux_fraction(
            p.d26Mg_river_water,
            d26Mg_silicate,
            p.d26Mg_carbonate,
        )
        flux_fraction_valid = np.isfinite(flux_fraction) and 0 <= flux_fraction <= 1

        if not flux_fraction_valid:
            flux_fraction = float("nan")
            f_silicate = float("nan")
            f_carbonate = float("nan")
            mass_balance = float("nan")
            if np.isfinite(f_mg):
                status = "incompatible_flux_endmembers"
                message = "River isotope value is outside the assumed source endmembers."
        else:
            f_silicate = p.F_river_total * flux_fraction
            f_carbonate = p.F_river_total - f_silicate
            reconstructed = (
                flux_fraction * d26Mg_silicate
                + (1.0 - flux_fraction) * p.d26Mg_carbonate
            )
            mass_balance = abs(p.d26Mg_river_water - reconstructed)

        within_calibration = status == "within_hu_calibration"
        mg_loss = 1.0 - f_mg if np.isfinite(f_mg) else float("nan")
        result = SilicateModelResult(
            success=True,
            message=message,
            data={
                "d26Mg_clay": float(d26Mg_clay),
                "f_Mg": f_mg,
                "f_Mg_unbounded": f_raw,
                "Mg_loss_fraction": mg_loss,
                "weathering_stage": stage,
                "model_status": status,
                "within_hu_calibration": within_calibration,
                "d26Mg_silicate": d26Mg_silicate,
                "d26Mg_carbonate": p.d26Mg_carbonate,
                "silicate_flux_fraction": flux_fraction,
                "F_silicate": f_silicate,
                "F_carbonate": f_carbonate,
                "F_river_total": p.F_river_total,
                "absolute_flux_is_conditional": True,
            },
            f_Mg=f_mg,
            Mg_loss_fraction=mg_loss,
            weathering_stage=stage,
            model_status=status,
            within_hu_calibration=within_calibration,
            d26Mg_silicate=d26Mg_silicate,
            d26Mg_carbonate=p.d26Mg_carbonate,
            silicate_flux_fraction=flux_fraction,
            F_silicate=f_silicate,
            F_carbonate=f_carbonate,
            SWI=flux_fraction * 100 if np.isfinite(flux_fraction) else float("nan"),
            mass_balance_check=mass_balance,
        )
        self._last_result = result
        return result

    def calculate_from_row(self, row_data: Dict[str, Any]) -> SilicateModelResult:
        possible_cols = ["delta_26_Mg_iso", "d26Mg", "δ26Mg", "delta_26Mg"]
        for col in possible_cols:
            if col in row_data and np.isfinite(float(row_data[col])):
                return self.calculate_weathering_flux(float(row_data[col]))
        return SilicateModelResult(
            success=False,
            message=f"No valid Mg isotope column found; supported names: {possible_cols}",
        )

    @staticmethod
    def _summary(samples: np.ndarray) -> Dict[str, Any]:
        samples = np.asarray(samples, dtype=float)
        samples = samples[np.isfinite(samples)]
        if samples.size == 0:
            return {
                "mean": np.nan,
                "std": np.nan,
                "median": np.nan,
                "ci_68": np.array([np.nan, np.nan]),
                "ci_95": np.array([np.nan, np.nan]),
                "samples": samples,
            }
        return {
            "mean": float(np.mean(samples)),
            "std": float(np.std(samples)),
            "median": float(np.median(samples)),
            "ci_68": np.percentile(samples, [16, 84]),
            "ci_95": np.percentile(samples, [2.5, 97.5]),
            "samples": samples,
        }

    def monte_carlo_analysis(
        self,
        d26Mg_clay: float,
        d26Mg_clay_std: Optional[float] = None,
        n_iterations: int = 10000,
        random_seed: Optional[int] = 42,
    ) -> Dict[str, Any]:
        """Propagate analytical and Hu-endmember uncertainties."""
        if n_iterations <= 0:
            raise ValueError("n_iterations must be positive")

        p = self.params
        clay_std = p.uncertainty_clay if d26Mg_clay_std is None else d26Mg_clay_std
        if clay_std < 0:
            raise ValueError("d26Mg_clay_std must be non-negative")

        rng = np.random.default_rng(random_seed)
        clay = rng.normal(d26Mg_clay, clay_std, n_iterations)
        protolith = self._draw_protolith(rng, n_iterations)
        delta = rng.uniform(*p.Delta_fluid_protolith_range, n_iterations)
        river = rng.normal(p.d26Mg_river_water, p.uncertainty_river, n_iterations)
        carbonate = rng.normal(
            p.d26Mg_carbonate,
            p.uncertainty_carbonate,
            n_iterations,
        )
        if p.uncertainty_river_flux_fraction > 0:
            river_flux = rng.normal(
                p.F_river_total,
                p.F_river_total * p.uncertainty_river_flux_fraction,
                n_iterations,
            )
        else:
            river_flux = np.full(n_iterations, p.F_river_total)

        alpha = 1.0 + delta / 1000.0
        initial_ratio = self._delta_to_relative_ratio(protolith)
        residue_ratio = self._delta_to_relative_ratio(clay)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            f_mg = (residue_ratio / initial_ratio) ** (1.0 / (alpha - 1.0))
        physical_f = np.isfinite(f_mg) & (f_mg >= 0) & (f_mg <= 1)

        released_ratio = np.full(n_iterations, np.nan)
        near_one = physical_f & np.isclose(f_mg, 1.0)
        released_ratio[near_one] = (
            alpha[near_one] * initial_ratio[near_one]
        )
        regular = physical_f & ~near_one
        released_ratio[regular] = initial_ratio[regular] * (
            (1.0 - f_mg[regular] ** alpha[regular])
            / (1.0 - f_mg[regular])
        )
        d26Mg_silicate = (released_ratio - 1.0) * 1000.0

        with np.errstate(invalid="ignore", divide="ignore"):
            flux_fraction = (river - carbonate) / (d26Mg_silicate - carbonate)
        valid_flux = (
            physical_f
            & np.isfinite(flux_fraction)
            & (flux_fraction >= 0)
            & (flux_fraction <= 1)
            & np.isfinite(river_flux)
            & (river_flux > 0)
        )

        f_samples = f_mg[physical_f]
        loss_samples = 1.0 - f_samples
        flux_fraction_samples = flux_fraction[valid_flux]
        f_silicate_samples = river_flux[valid_flux] * flux_fraction_samples
        f_carbonate_samples = river_flux[valid_flux] - f_silicate_samples

        return {
            "d26Mg_protolith": self._summary(protolith),
            "f_Mg": self._summary(f_samples),
            "Mg_loss_fraction": self._summary(loss_samples),
            "d26Mg_silicate": self._summary(d26Mg_silicate[valid_flux]),
            "silicate_flux_fraction": self._summary(flux_fraction_samples),
            "F_silicate": self._summary(f_silicate_samples),
            "F_carbonate": self._summary(f_carbonate_samples),
            # Backward-compatible percentage view of dissolved Mg partition.
            "SWI": self._summary(flux_fraction_samples * 100.0),
            "n_physical_weathering": int(np.sum(physical_f)),
            "n_valid_flux": int(np.sum(valid_flux)),
            "n_total": int(n_iterations),
            "valid_weathering_fraction": float(np.mean(physical_f)),
            "valid_flux_fraction": float(np.mean(valid_flux)),
            "protolith_prior_range": p.d26Mg_UCC_range,
        }

    def monte_carlo_profile(
        self,
        d26Mg_clay: np.ndarray,
        d26Mg_clay_std: Optional[np.ndarray] = None,
        baseline_count: int = 5,
        baseline_position: str = "start",
        n_iterations: int = 10000,
        random_seed: Optional[int] = 42,
    ) -> Dict[str, Any]:
        """Calculate profile weathering proxies with shared nuisance parameters.

        The primary profile result is the baseline-free Mg release fraction,
        ``1 - f_Mg``. Under constant rock supply and protolith Mg content this
        is proportional to weathering Mg flux. A two-state L1 change point is
        detected automatically, and the deep/shallow state ratio quantifies
        the fold change without a manually selected reference interval.

        Baseline-normalized fields are retained for API compatibility, but are
        not required by the release-fraction or change-point calculations.
        """
        values = np.asarray(d26Mg_clay, dtype=float)
        if values.ndim != 1 or values.size == 0:
            raise ValueError("d26Mg_clay must be a non-empty one-dimensional array")
        if not 1 <= baseline_count <= values.size:
            raise ValueError("baseline_count must be between 1 and the sample count")
        if baseline_position not in ("start", "end"):
            raise ValueError("baseline_position must be 'start' or 'end'")
        if n_iterations <= 0:
            raise ValueError("n_iterations must be positive")

        p = self.params
        if d26Mg_clay_std is None:
            std = np.full(values.size, p.uncertainty_clay)
        else:
            std = np.asarray(d26Mg_clay_std, dtype=float)
            if std.shape != values.shape:
                raise ValueError("d26Mg_clay_std must match d26Mg_clay")
            std = np.where(np.isfinite(std), std, p.uncertainty_clay)
        if np.any(std < 0):
            raise ValueError("d26Mg_clay_std must be non-negative")

        rng = np.random.default_rng(random_seed)
        clay = rng.normal(values[:, None], std[:, None], (values.size, n_iterations))
        protolith = self._draw_protolith(rng, n_iterations)
        delta = rng.uniform(*p.Delta_fluid_protolith_range, n_iterations)
        river = rng.normal(p.d26Mg_river_water, p.uncertainty_river, n_iterations)
        carbonate = rng.normal(
            p.d26Mg_carbonate,
            p.uncertainty_carbonate,
            n_iterations,
        )
        if p.uncertainty_river_flux_fraction > 0:
            river_flux = rng.normal(
                p.F_river_total,
                p.F_river_total * p.uncertainty_river_flux_fraction,
                n_iterations,
            )
        else:
            river_flux = np.full(n_iterations, p.F_river_total)

        alpha = 1.0 + delta / 1000.0
        initial_ratio = self._delta_to_relative_ratio(protolith)
        residue_ratio = self._delta_to_relative_ratio(clay)
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            f_mg = (residue_ratio / initial_ratio[None, :]) ** (
                1.0 / (alpha[None, :] - 1.0)
            )
        physical_f = np.isfinite(f_mg) & (f_mg >= 0) & (f_mg <= 1)

        released_ratio = np.full_like(f_mg, np.nan)
        near_one = physical_f & np.isclose(f_mg, 1.0)
        alpha_grid = np.broadcast_to(alpha, f_mg.shape)
        initial_grid = np.broadcast_to(initial_ratio, f_mg.shape)
        released_ratio[near_one] = (
            alpha_grid[near_one] * initial_grid[near_one]
        )
        regular = physical_f & ~near_one
        released_ratio[regular] = initial_grid[regular] * (
            (1.0 - f_mg[regular] ** alpha_grid[regular])
            / (1.0 - f_mg[regular])
        )
        d26Mg_silicate = (released_ratio - 1.0) * 1000.0

        with np.errstate(invalid="ignore", divide="ignore"):
            flux_fraction = (
                (river[None, :] - carbonate[None, :])
                / (d26Mg_silicate - carbonate[None, :])
            )
        valid_flux = (
            physical_f
            & np.isfinite(flux_fraction)
            & (flux_fraction >= 0)
            & (flux_fraction <= 1)
            & np.isfinite(river_flux[None, :])
            & (river_flux[None, :] > 0)
        )
        f_silicate = np.where(
            valid_flux,
            flux_fraction * river_flux[None, :],
            np.nan,
        )

        baseline_flux = np.full(n_iterations, np.nan)
        if baseline_position == "start":
            baseline_block = f_silicate[:baseline_count, :]
            baseline_indices = np.arange(baseline_count)
        else:
            baseline_block = f_silicate[-baseline_count:, :]
            baseline_indices = np.arange(values.size - baseline_count, values.size)
        for draw in range(n_iterations):
            finite = baseline_block[:, draw][np.isfinite(baseline_block[:, draw])]
            if finite.size:
                baseline_flux[draw] = np.median(finite)

        with np.errstate(invalid="ignore", divide="ignore"):
            multiplier = f_silicate / baseline_flux[None, :]
        multiplier[~np.isfinite(multiplier)] = np.nan

        # Under constant rock supply and protolith Mg concentration, Mg loss
        # per unit rock is proportional to the total weathering Mg flux.
        mg_loss = np.where(physical_f, 1.0 - f_mg, np.nan)
        baseline_loss = np.full(n_iterations, np.nan)
        baseline_loss_block = mg_loss[baseline_indices, :]
        for draw in range(n_iterations):
            finite = baseline_loss_block[:, draw][
                np.isfinite(baseline_loss_block[:, draw])
            ]
            if finite.size:
                baseline_loss[draw] = np.median(finite)
        with np.errstate(invalid="ignore", divide="ignore"):
            weathering_multiplier = mg_loss / baseline_loss[None, :]
        weathering_multiplier[
            ~np.isfinite(weathering_multiplier) | (weathering_multiplier < 0)
        ] = np.nan

        medians = np.full(values.size, np.nan)
        ci_low = np.full(values.size, np.nan)
        ci_high = np.full(values.size, np.nan)
        valid_fraction = np.zeros(values.size)
        weathering_medians = np.full(values.size, np.nan)
        weathering_ci_low = np.full(values.size, np.nan)
        weathering_ci_high = np.full(values.size, np.nan)
        weathering_valid_fraction = np.zeros(values.size)
        release_medians = np.full(values.size, np.nan)
        release_ci_low = np.full(values.size, np.nan)
        release_ci_high = np.full(values.size, np.nan)
        release_valid_fraction = np.zeros(values.size)
        conditional_flux_medians = np.full(values.size, np.nan)
        conditional_flux_ci_low = np.full(values.size, np.nan)
        conditional_flux_ci_high = np.full(values.size, np.nan)
        conditional_flux_valid_fraction = np.zeros(values.size)
        for index in range(values.size):
            summary = self._summary(multiplier[index])
            medians[index] = summary["median"]
            ci_low[index], ci_high[index] = summary["ci_95"]
            valid_fraction[index] = np.mean(np.isfinite(multiplier[index]))
            weathering_summary = self._summary(weathering_multiplier[index])
            weathering_medians[index] = weathering_summary["median"]
            weathering_ci_low[index], weathering_ci_high[index] = (
                weathering_summary["ci_95"]
            )
            weathering_valid_fraction[index] = np.mean(
                np.isfinite(weathering_multiplier[index])
            )
            release_summary = self._summary(mg_loss[index])
            release_medians[index] = release_summary["median"]
            release_ci_low[index], release_ci_high[index] = release_summary["ci_95"]
            release_valid_fraction[index] = np.mean(np.isfinite(mg_loss[index]))
            flux_summary = self._summary(f_silicate[index])
            conditional_flux_medians[index] = flux_summary["median"]
            conditional_flux_ci_low[index], conditional_flux_ci_high[index] = (
                flux_summary["ci_95"]
            )
            conditional_flux_valid_fraction[index] = np.mean(
                np.isfinite(f_silicate[index])
            )

        def _find_l1_split(profile):
            finite_profile = profile[np.isfinite(profile)]
            if finite_profile.size < 4:
                return None, np.nan, np.nan, np.nan, np.nan
            single_cost = np.sum(
                np.abs(finite_profile - np.median(finite_profile))
            )
            best = None
            for split in range(2, values.size - 1):
                shallow = profile[:split]
                shallow = shallow[np.isfinite(shallow)]
                deep = profile[split:]
                deep = deep[np.isfinite(deep)]
                if shallow.size < 2 or deep.size < 2:
                    continue
                shallow_median = np.median(shallow)
                deep_median = np.median(deep)
                cost = (
                    np.sum(np.abs(shallow - shallow_median))
                    + np.sum(np.abs(deep - deep_median))
                )
                if best is None or cost < best[0]:
                    best = (cost, split, shallow_median, deep_median)
            if best is None:
                return None, np.nan, np.nan, np.nan, np.nan
            cost, split, shallow_median, deep_median = best
            improvement = (
                1.0 - cost / single_cost if single_cost > 0 else np.nan
            )
            return split, shallow_median, deep_median, cost, improvement

        # Use the complete posterior-median profile for the reported split.
        # This prevents draws with several boundary-incompatible samples from
        # moving the primary split merely because their values are missing.
        best_split, _, _, _, central_improvement = _find_l1_split(
            release_medians
        )

        split_points = []
        cost_improvements = []
        for draw in range(n_iterations):
            profile = mg_loss[:, draw]
            split, _, _, _, improvement = _find_l1_split(profile)
            if split is None:
                continue
            split_points.append(split)
            cost_improvements.append(improvement)

        shallow_states = []
        deep_states = []
        deep_to_shallow_ratios = []
        if best_split is not None:
            for draw in range(n_iterations):
                shallow = mg_loss[:best_split, draw]
                shallow = shallow[np.isfinite(shallow)]
                deep = mg_loss[best_split:, draw]
                deep = deep[np.isfinite(deep)]
                if shallow.size < 2 or deep.size < 2:
                    continue
                shallow_median = np.median(shallow)
                deep_median = np.median(deep)
                shallow_states.append(shallow_median)
                deep_states.append(deep_median)
                deep_to_shallow_ratios.append(
                    deep_median / shallow_median
                    if shallow_median > 0
                    else np.nan
                )

        split_points = np.asarray(split_points, dtype=int)
        if split_points.size:
            split_mode = int(np.bincount(split_points).argmax())
        else:
            split_mode = None
        split_summary = self._summary(split_points.astype(float))

        return {
            "d26Mg_protolith": self._summary(protolith),
            "flux_multiplier_median": medians,
            "flux_multiplier_ci95_low": ci_low,
            "flux_multiplier_ci95_high": ci_high,
            "valid_multiplier_fraction": valid_fraction,
            "weathering_flux_multiplier_median": weathering_medians,
            "weathering_flux_multiplier_ci95_low": weathering_ci_low,
            "weathering_flux_multiplier_ci95_high": weathering_ci_high,
            "valid_weathering_multiplier_fraction": weathering_valid_fraction,
            "Mg_release_fraction_median": release_medians,
            "Mg_release_fraction_ci95_low": release_ci_low,
            "Mg_release_fraction_ci95_high": release_ci_high,
            "valid_Mg_release_fraction": release_valid_fraction,
            "conditional_F_silicate_median": conditional_flux_medians,
            "conditional_F_silicate_ci95_low": conditional_flux_ci_low,
            "conditional_F_silicate_ci95_high": conditional_flux_ci_high,
            "valid_conditional_F_silicate": conditional_flux_valid_fraction,
            "baseline_flux": self._summary(baseline_flux),
            "baseline_Mg_loss_fraction": self._summary(baseline_loss),
            "baseline_count": int(baseline_count),
            "baseline_position": baseline_position,
            "baseline_indices": baseline_indices,
            "n_iterations": int(n_iterations),
            "absolute_flux_scale_cancels": True,
            "weathering_flux_assumes_constant_rock_supply_and_Mg": True,
            "protolith_is_shared_across_profile": True,
            "protolith_prior_range": p.d26Mg_UCC_range,
            "change_point": {
                "after_sample_best": best_split,
                "after_sample_mode": split_mode,
                "after_sample_median": split_summary["median"],
                "after_sample_ci95": split_summary["ci_95"],
                "shallow_state": self._summary(np.asarray(shallow_states)),
                "deep_state": self._summary(np.asarray(deep_states)),
                "deep_to_shallow_ratio": self._summary(
                    np.asarray(deep_to_shallow_ratios)
                ),
                "cost_improvement": self._summary(np.asarray(cost_improvements)),
                "central_cost_improvement": central_improvement,
                "probability_deep_greater": float(
                    np.mean(np.asarray(deep_states) > np.asarray(shallow_states))
                ) if deep_states else np.nan,
                "n_valid": int(split_points.size),
                "n_valid_state_ratios": int(
                    np.sum(np.isfinite(deep_to_shallow_ratios))
                ),
            },
        }

    def simulate_transient(
        self,
        d26Mg_clay_timeseries: np.ndarray,
        times: np.ndarray,
    ) -> list:
        if len(d26Mg_clay_timeseries) != len(times):
            raise ValueError("d26Mg_clay_timeseries and times must have equal length")
        return [self.calculate_weathering_flux(value) for value in d26Mg_clay_timeseries]

    def get_info(self) -> Dict[str, Any]:
        p = self.params
        return {
            "component_type": self.COMPONENT_TYPE,
            "basin": self.basin,
            "reference": "Hu et al. (2023)",
            "description": "Siliciclastic Mg-isotope weathering and conditional flux model",
            "end_members": {
                "protolith": {
                    "d26Mg_reference": p.d26Mg_UCC,
                    "d26Mg_uniform_prior": p.d26Mg_UCC_range,
                    "shared_across_profile": True,
                },
                "carbonate_flux": {"d26Mg": p.d26Mg_carbonate},
                "river_water": {"d26Mg": p.d26Mg_river_water},
            },
            "fractionation_factors": {
                "Delta_fluid_protolith": p.Delta_fluid_protolith,
                "range": p.Delta_fluid_protolith_range,
            },
            "flux_prior": {
                "F_river_total": p.F_river_total,
                "absolute_flux_is_conditional": True,
            },
            "applicability": [
                "clay-sized terrigenous siliciclastic sediment",
                "incipient-to-intermediate silicate weathering",
                "advanced weathering only as an explicitly marked extrapolation",
            ],
        }


SilicateWeatheringModel = SilicateMgSystem
