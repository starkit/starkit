from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
from starkit.base.operations.base import SpectralOperationModel
from starkit.base.operations.spectrograph import (
    InstrumentConvolveGrating,
    Interpolate,
    Normalize,
    NormalizeParts,
    SpectrographOperationModel,
)
from starkit.base.operations.stellar import RotationalBroadening
from starkit.base.operations.stellar import DopplerShift
from starkit.fitkit import priors
from starkit.fitkit.likelihoods import SpectralChi2Likelihood, SpectralL1Likelihood
from starkit.gridkit import load_grid
from starkit.utils.spectral import prepare_observed

from .config import AnalysisConfig, PriorConfig
from .mpi import RuntimeContext

LOGGER = logging.getLogger(__name__)


@dataclass(frozen=True)
class BuiltModel:
    grid: object
    model: object
    fit_model: object
    prior_objects: list[object]
    prior_names: list[str]
    parameter_name_map: dict[str, str]


class CropToObservedRange(SpectrographOperationModel):
    requires_observed_spectrum = True
    operation_name = "crop_to_observed_range"

    def __init__(self, observed) -> None:
        super(SpectralOperationModel, self).__init__()
        self.observed = prepare_observed(observed)

    def evaluate(self, wavelength, flux):
        lower = self.observed.wavelength.value.min()
        upper = self.observed.wavelength.value.max()
        mask = (wavelength >= lower) & (wavelength <= upper)
        return wavelength[mask], flux[mask]


def load_analysis_grid(config: AnalysisConfig):
    grid = load_grid(
        str(config.run.grid_path),
        wavelength_type=config.run.grid_wavelength_type,
    )
    LOGGER.info(
        "Loaded grid %s with R=%s R_sampling=%s wavelength_type=%s native_wavelength_type=%s",
        config.run.grid_path,
        getattr(grid, "R", None),
        getattr(grid, "R_sampling", None),
        getattr(grid, "wavelength_type", None),
        grid.meta_grid.get("wavelength_type"),
    )
    return grid


def build_model(config: AnalysisConfig, grid, observed, *, normalization_parts: list[np.ndarray] | None = None) -> BuiltModel:
    _validate_resolution_configuration(config, grid)
    model = grid

    if config.model.rotational_broadening:
        model = model | RotationalBroadening.from_grid(
            grid,
            vrot=np.array([config.model.rotational_velocity]),
        )
    if config.model.doppler_shift:
        model = model | DopplerShift(vrad=config.model.radial_velocity)
    if config.model.instrument_resolution is not None:
        model = model | InstrumentConvolveGrating.from_grid(
            grid,
            R=config.model.instrument_resolution,
        )
    if config.model.crop_to_observed_range:
        model = model | CropToObservedRange(observed)
    if config.model.interpolate:
        model = model | Interpolate(observed)
    normalizer = _build_normalizer(config, observed, normalization_parts)
    if normalizer is not None:
        model = model | normalizer

    likelihood = _build_likelihood(config.model.likelihood, observed)
    fit_model = model | likelihood
    parameter_name_map = _build_parameter_name_map(fit_model.param_names)
    prior_objects, prior_names, resolved_prior_names = _configure_fit_priors(
        fit_model,
        config.priors,
        parameter_name_map,
    )
    LOGGER.info(
        "Configured priors: %s",
        ", ".join(f"{user}->{actual}" for user, actual in zip(prior_names, resolved_prior_names)),
    )
    return BuiltModel(
        grid=grid,
        model=model,
        fit_model=fit_model,
        prior_objects=prior_objects,
        prior_names=prior_names,
        parameter_name_map=parameter_name_map,
    )


def run_fit(built_model: BuiltModel, runtime: RuntimeContext, clean_up: bool):
    from starkit.fitkit.multinest.base import MultiNest

    fit = MultiNest(
        built_model.fit_model,
        built_model.prior_objects,
        run_dir=str(runtime.run_dir),
    )
    return fit.run(clean_up=(clean_up and runtime.is_root))


def evaluate_best_fit(model, result):
    for name, value in result.median.items():
        setattr(model, name, value)
    return model()


def _build_likelihood(name: str, observed):
    normalized = name.lower()
    if normalized == "chi2":
        return SpectralChi2Likelihood(observed)
    if normalized == "l1":
        return SpectralL1Likelihood(observed)
    raise ValueError(f"Unsupported likelihood: {name}")


def _build_prior(config: PriorConfig):
    kind = config.type.lower()
    if kind == "fixed":
        if config.value is None:
            raise ValueError("Fixed prior requires 'value'")
        return priors.FixedPrior(config.value)
    if kind == "uniform":
        if config.min is None or config.max is None:
            raise ValueError("Uniform prior requires 'min' and 'max'")
        return priors.UniformPrior(config.min, config.max)
    if kind == "gaussian":
        if config.mean is None or config.sigma is None:
            raise ValueError("Gaussian prior requires 'mean' and 'sigma'")
        return priors.GaussianPrior(config.mean, config.sigma)
    raise ValueError(f"Unsupported prior type: {config.type}")


def _configure_fit_priors(
    fit_model,
    configured_priors: dict[str, PriorConfig],
    parameter_name_map: dict[str, str],
) -> tuple[list[object], list[str], list[str]]:
    resolved_lookup: dict[str, tuple[str, object]] = {}
    for user_name, prior_config in configured_priors.items():
        resolved_name = parameter_name_map.get(user_name)
        if resolved_name is None:
            raise ValueError(
                f"Configured prior '{user_name}' is not a parameter on the fit model. "
                f"Available parameters: {', '.join(sorted(parameter_name_map))}"
            )
        if resolved_name in resolved_lookup:
            previous_user_name, _ = resolved_lookup[resolved_name]
            raise ValueError(
                f"Configured priors '{previous_user_name}' and '{user_name}' resolve to the same "
                f"fit-model parameter '{resolved_name}'."
            )
        resolved_lookup[resolved_name] = (user_name, _build_prior(prior_config))

    for param_name in fit_model.param_names:
        fit_model.fixed[param_name] = True
    for resolved_name in resolved_lookup:
        fit_model.fixed[resolved_name] = False

    ordered_prior_objects: list[object] = []
    ordered_user_names: list[str] = []
    ordered_actual_names: list[str] = []
    for param_name in fit_model.param_names:
        if fit_model.fixed[param_name]:
            continue
        user_name, prior_object = resolved_lookup[param_name]
        ordered_prior_objects.append(prior_object)
        ordered_user_names.append(user_name)
        ordered_actual_names.append(param_name)

    return ordered_prior_objects, ordered_user_names, ordered_actual_names


def get_grid_sampling(grid) -> tuple[float | None, float | None]:
    return getattr(grid, "R", None), getattr(grid, "R_sampling", None)


def _build_normalizer(
    config: AnalysisConfig,
    observed,
    normalization_parts: list[np.ndarray] | None = None,
):
    npol = config.model.normalize_npol
    if npol is None:
        return None
    if config.model.normalize_parts is None:
        return Normalize(observed, npol)
    if config.model.normalize_parts in {"segments", "autodetect"}:
        if normalization_parts is None:
            raise ValueError(f"model.normalize_parts={config.model.normalize_parts!r} requires computed normalization parts.")
        return NormalizeParts(observed, normalization_parts, npol)
    raise ValueError(f"Unsupported model.normalize_parts: {config.model.normalize_parts!r}")


def _build_parameter_name_map(param_names) -> dict[str, str]:
    mapping: dict[str, str] = {}
    for param_name in param_names:
        base_name = param_name.rsplit("_", 1)[0]
        if base_name in mapping and mapping[base_name] != param_name:
            raise ValueError(
                f"Ambiguous compound-model parameter base name '{base_name}' "
                f"for parameters '{mapping[base_name]}' and '{param_name}'"
            )
        mapping[base_name] = param_name
    return mapping


def _validate_resolution_configuration(config: AnalysisConfig, grid) -> None:
    grid_resolution = getattr(grid, "R", None)
    requested_resolution = config.model.instrument_resolution

    if requested_resolution is not None and grid_resolution is not None:
        if float(requested_resolution) >= float(grid_resolution):
            raise ValueError(
                "model.instrument_resolution must be smaller than the grid resolution "
                f"(requested {requested_resolution}, grid {grid_resolution})"
            )

    resolution_prior = config.priors.get("R")
    if resolution_prior is None or grid_resolution is None:
        return

    if resolution_prior.type.lower() == "uniform" and resolution_prior.max is not None:
        if float(resolution_prior.max) >= float(grid_resolution):
            raise ValueError(
                "Uniform prior for R must have max smaller than the grid resolution "
                f"(prior max {resolution_prior.max}, grid {grid_resolution})"
            )
    if resolution_prior.type.lower() == "fixed" and resolution_prior.value is not None:
        if float(resolution_prior.value) >= float(grid_resolution):
            raise ValueError(
                "Fixed prior for R must be smaller than the grid resolution "
                f"(prior value {resolution_prior.value}, grid {grid_resolution})"
            )
