from __future__ import annotations

import importlib
import json
import logging
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter

from .config import AnalysisConfig
from .spectrum import convert_wavelength_values
from starkit.evolutionkit.yyisochrones_galactic_center import predict_logg
from starkit.evolutionkit.yyisochrones_galactic_center import predict_logg_per_slice

LOGGER = logging.getLogger(__name__)
FIXED_PARAMETER_TOLERANCE = 1e-12
CORNER_MIN_SAMPLES_FOR_CONTOURS = 32
CORNER_LABELS = {
    "R": "R",
    "afe": "[alpha/Fe]",
    "alpha": "[alpha/Fe]",
    "alphafe": "[alpha/Fe]",
    "cfe": "[C/Fe]",
    "logg": "log g",
    "mh": "[M/H]",
    "nfe": "[N/Fe]",
    "ofe": "[O/Fe]",
    "ratio": "C12/C13",
    "teff": "Teff",
    "vrad": "vrad",
    "vrot": "vrot",
}


def write_outputs(
    config: AnalysisConfig,
    result,
    observed,
    best_fit_wavelength,
    best_fit_flux,
    built_model=None,
    include_regions: np.ndarray | None = None,
    exclude_regions: np.ndarray | None = None,
) -> None:
    output_dir = config.outputs.directory
    output_dir.mkdir(parents=True, exist_ok=True)
    source_wavelength_type, output_wavelength_type = _resolve_output_wavelength_types(config, built_model)
    observed_wavelength = _convert_output_wavelength(
        observed.wavelength.value,
        config,
        source_wavelength_type,
        output_wavelength_type,
    )
    best_fit_output_wavelength = _convert_output_wavelength(
        best_fit_wavelength,
        config,
        source_wavelength_type,
        output_wavelength_type,
    )
    output_include_regions = _convert_output_regions(
        include_regions,
        config,
        source_wavelength_type,
        output_wavelength_type,
    )
    output_exclude_regions = _convert_output_regions(
        exclude_regions,
        config,
        source_wavelength_type,
        output_wavelength_type,
    )

    if config.outputs.write_hdf5:
        hdf5_path = output_dir / _filename(config, "hdf5", f"{config.run.name}_result.h5")
        result.to_hdf(str(hdf5_path))
        LOGGER.info("Wrote %s", hdf5_path)

    if config.outputs.write_bestfit_spectrum:
        best_fit_path = output_dir / _filename(
            config,
            "bestfit_spectrum",
            f"{config.run.name}_bestfit_spectrum.txt",
        )
        best_fit_data = np.column_stack([best_fit_output_wavelength, best_fit_flux])
        np.savetxt(best_fit_path, best_fit_data, header="wavelength flux")
        LOGGER.info("Wrote %s", best_fit_path)

    if config.outputs.write_masked_observed:
        observed_path = output_dir / _filename(
            config,
            "masked_observed",
            f"{config.run.name}_masked_observed.txt",
        )
        observed_uncertainty = getattr(observed, "uncertainty", None)
        if observed_uncertainty is None:
            observed_data = np.column_stack([observed_wavelength, observed.flux.value])
            header = "wavelength flux"
        else:
            observed_data = np.column_stack(
                [observed_wavelength, observed.flux.value, observed_uncertainty.value]
            )
            header = "wavelength flux uncertainty"
        np.savetxt(observed_path, observed_data, header=header)
        LOGGER.info("Wrote %s", observed_path)

    if config.outputs.write_run_summary:
        summary_path = output_dir / _filename(
            config,
            "run_summary",
            f"{config.run.name}_summary.json",
        )
        grid_parameter_order = _get_grid_parameter_order(getattr(built_model, "grid", None))
        parameter_name_map = getattr(built_model, "parameter_name_map", None)
        median_summary = _order_summary_parameters(
            {key: _json_value(value) for key, value in result.median.items()},
            grid_parameter_order=grid_parameter_order,
            parameter_name_map=parameter_name_map,
        )
        mean_summary = _order_summary_parameters(
            {key: _json_value(value) for key, value in result.mean.items()},
            grid_parameter_order=grid_parameter_order,
            parameter_name_map=parameter_name_map,
        )
        derived = _build_derived_summary(
            median=median_summary,
            mean=mean_summary,
            parameter_name_map=parameter_name_map,
        )
        summary = {
            "median": median_summary,
            "mean": mean_summary,
        }
        if derived:
            summary["derived"] = derived
        with summary_path.open("w", encoding="ascii") as handle:
            json.dump(summary, handle, indent=2)
            handle.write("\n")
        LOGGER.info("Wrote %s", summary_path)

    if config.outputs.write_plot:
        plot_path = output_dir / _filename(config, "plot", f"{config.run.name}_fit.pdf")
        figure, axis = plt.subplots(figsize=(14, 7))
        plot_min, plot_max = _resolve_plot_bounds(observed_wavelength, best_fit_output_wavelength)
        _add_region_spans(
            axis,
            output_include_regions,
            plot_min,
            plot_max,
            color="#2f9e44",
            alpha=0.12,
            label="included region",
        )
        _add_region_spans(
            axis,
            output_exclude_regions,
            plot_min,
            plot_max,
            color="#c92a2a",
            alpha=0.22,
            label="excluded region",
        )
        axis.plot(observed_wavelength, observed.flux.value, label="observed")
        axis.plot(best_fit_output_wavelength, best_fit_flux, alpha=0.7, label="best fit")
        axis.set_xlim(plot_min, plot_max)
        axis.set_xlabel(_wavelength_axis_label(observed.wavelength.unit, output_wavelength_type))
        axis.set_ylabel(f"Flux ({observed.flux.unit})")
        axis.legend(loc="best")
        figure.tight_layout()
        figure.savefig(plot_path, dpi=200)
        plt.close(figure)
        LOGGER.info("Wrote %s", plot_path)

    if config.outputs.write_corner_plot:
        corner_path = output_dir / _filename(
            config,
            "corner_plot",
            f"{config.run.name}_corner.pdf",
        )
        if _write_corner_plot(
            corner_path,
            result=result,
            built_model=built_model,
        ):
            LOGGER.info("Wrote %s", corner_path)


def _filename(config: AnalysisConfig, key: str, default: str) -> str:
    return config.outputs.filenames.get(key, default)


def _resolve_output_wavelength_types(
    config: AnalysisConfig,
    built_model=None,
) -> tuple[str | None, str | None]:
    source_wavelength_type = _resolve_analysis_wavelength_type(config, built_model)
    output_wavelength_type = source_wavelength_type
    if config.outputs.wavelength_type is not None:
        output_wavelength_type = config.outputs.wavelength_type
    if config.outputs.wavelength_type is not None and source_wavelength_type is None:
        raise ValueError(
            "outputs.wavelength_type requires the analysis wavelength type to be known. "
            "Set run.grid_wavelength_type or input_spectrum.wavelength_type."
        )
    return source_wavelength_type, output_wavelength_type


def _resolve_analysis_wavelength_type(config: AnalysisConfig, built_model=None) -> str | None:
    grid = getattr(built_model, "grid", None)
    grid_wavelength_type = getattr(grid, "wavelength_type", None)
    if grid_wavelength_type is not None:
        return grid_wavelength_type
    if config.run.grid_wavelength_type is not None:
        return config.run.grid_wavelength_type
    return config.input_spectrum.wavelength_type


def _convert_output_wavelength(
    values: np.ndarray,
    config: AnalysisConfig,
    source_wavelength_type: str | None,
    output_wavelength_type: str | None,
) -> np.ndarray:
    return convert_wavelength_values(
        np.asarray(values, dtype=float),
        config.input_spectrum.wavelength_unit,
        source_wavelength_type,
        output_wavelength_type,
    )


def _convert_output_regions(
    regions: np.ndarray | None,
    config: AnalysisConfig,
    source_wavelength_type: str | None,
    output_wavelength_type: str | None,
) -> np.ndarray | None:
    if regions is None:
        return None
    converted = _convert_output_wavelength(regions, config, source_wavelength_type, output_wavelength_type)
    region_array = np.asarray(converted, dtype=float).reshape(-1, 2)
    return np.asarray(
        [[min(start, stop), max(start, stop)] for start, stop in region_array],
        dtype=float,
    ).reshape(-1, 2)


def _wavelength_axis_label(unit, wavelength_type: str | None) -> str:
    if wavelength_type is None:
        return f"Wavelength ({unit})"
    return f"Wavelength ({unit}, {wavelength_type})"


def _json_value(value):
    if isinstance(value, np.generic):
        return value.item()
    return value


def _resolve_plot_bounds(observed_wavelength: np.ndarray, best_fit_wavelength: np.ndarray) -> tuple[float, float]:
    wavelength_arrays = [
        np.asarray(observed_wavelength, dtype=float).ravel(),
        np.asarray(best_fit_wavelength, dtype=float).ravel(),
    ]
    finite_arrays = [values[np.isfinite(values)] for values in wavelength_arrays if values.size > 0]
    finite_arrays = [values for values in finite_arrays if values.size > 0]
    if not finite_arrays:
        raise ValueError("Cannot determine spectrum plot bounds from empty wavelength arrays")

    combined = np.concatenate(finite_arrays)
    return float(np.min(combined)), float(np.max(combined))


def _add_region_spans(
    axis,
    regions: np.ndarray | None,
    plot_min: float,
    plot_max: float,
    *,
    color: str,
    alpha: float,
    label: str,
) -> None:
    if regions is None:
        return

    region_array = np.asarray(regions, dtype=float).reshape(-1, 2)
    label_used = False
    for start, stop in region_array:
        clipped_start = max(float(start), plot_min)
        clipped_stop = min(float(stop), plot_max)
        if clipped_start > clipped_stop:
            continue
        if clipped_start == clipped_stop:
            axis.axvline(
                clipped_start,
                color=color,
                alpha=alpha,
                zorder=0,
                label=label if not label_used else None,
            )
            label_used = True
            continue
        axis.axvspan(
            clipped_start,
            clipped_stop,
            color=color,
            alpha=alpha,
            zorder=0,
            label=label if not label_used else None,
        )
        label_used = True


def _write_corner_plot(
    path: Path,
    *,
    result,
    built_model=None,
) -> bool:
    posterior_data = getattr(result, "posterior_data", None)
    if posterior_data is None:
        LOGGER.info("Skipping corner plot because the fit result does not expose posterior samples.")
        return False

    parameter_names = _resolve_corner_parameter_names(result, posterior_data, built_model)
    if not parameter_names:
        LOGGER.info("Skipping corner plot because no posterior parameters were available.")
        return False

    samples = posterior_data.loc[:, parameter_names].to_numpy(dtype=float, copy=True)
    if samples.size == 0:
        LOGGER.info("Skipping corner plot because the posterior sample table is empty.")
        return False

    weights = None
    if "weights" in posterior_data.columns:
        weights = np.asarray(posterior_data["weights"], dtype=float)

    finite_mask = np.all(np.isfinite(samples), axis=1)
    if weights is not None:
        finite_mask &= np.isfinite(weights)

    samples = samples[finite_mask]
    if weights is not None:
        weights = weights[finite_mask]

    if samples.size == 0:
        LOGGER.warning("Skipping corner plot because no finite posterior samples remain after filtering.")
        return False

    varying_mask = np.ptp(samples, axis=0) > FIXED_PARAMETER_TOLERANCE
    if not varying_mask.any():
        LOGGER.warning("Skipping corner plot because all posterior parameters are effectively fixed.")
        return False

    dropped = [name for name, keep in zip(parameter_names, varying_mask) if not keep]
    if dropped:
        LOGGER.info("Dropping fixed parameter(s) from corner plot: %s", ", ".join(dropped))

    plot_parameter_names = [name for name, keep in zip(parameter_names, varying_mask) if keep]
    plot_samples = samples[:, varying_mask]
    plot_labels = [
        _format_corner_label(name, getattr(built_model, "parameter_name_map", None))
        for name in plot_parameter_names
    ]

    if weights is not None:
        positive_weight_mask = weights > 0.0
        if not positive_weight_mask.any():
            LOGGER.warning("Skipping corner plot because posterior weights are not positive.")
            return False
        plot_samples = plot_samples[positive_weight_mask]
        weights = weights[positive_weight_mask]

    plot_contours = _should_plot_corner_contours(plot_samples, weights)

    corner_module = _load_corner_module()
    if corner_module is None:
        LOGGER.warning("Skipping corner plot because the 'corner' package is not installed.")
        return False

    figure_size = max(6.0, 2.8 * len(plot_parameter_names))
    figure = plt.figure(figsize=(figure_size, figure_size))
    title_kwargs = {"fontsize": 10}
    figure = corner_module.corner(
        plot_samples,
        labels=plot_labels,
        weights=weights,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_fmt=".3g",
        title_kwargs=title_kwargs,
        label_kwargs={"fontsize": 12},
        plot_contours=plot_contours,
        fill_contours=plot_contours,
        quiet=True,
        fig=figure,
    )
    _apply_corner_teff_formatting(
        figure,
        parameter_names=plot_parameter_names,
        labels=plot_labels,
        samples=plot_samples,
        weights=weights,
        parameter_name_map=getattr(built_model, "parameter_name_map", None),
        title_kwargs=title_kwargs,
    )
    figure.suptitle("Posterior Distributions", fontsize=16)
    figure.subplots_adjust(top=0.95)
    for axis in figure.axes:
        axis.tick_params(axis="both", labelsize=8)
    figure.savefig(path, dpi=200)
    plt.close(figure)
    return True


def _resolve_corner_parameter_names(result, posterior_data, built_model=None) -> list[str]:
    parameter_names = _resolve_posterior_parameter_names(result, posterior_data)
    if not parameter_names:
        return []

    grid_parameter_order = _get_grid_parameter_order(getattr(built_model, "grid", None))
    parameter_name_map = getattr(built_model, "parameter_name_map", None)
    if not grid_parameter_order:
        return parameter_names

    selected: list[str] = []
    resolved_parameter_name_map = parameter_name_map or {}
    for grid_parameter in grid_parameter_order:
        actual_name = resolved_parameter_name_map.get(grid_parameter, grid_parameter)
        if actual_name in parameter_names and actual_name not in selected:
            selected.append(actual_name)

    if selected:
        return selected

    LOGGER.info("Falling back to all posterior parameters because none matched the grid parameter list.")
    return parameter_names


def _resolve_posterior_parameter_names(result, posterior_data) -> list[str]:
    parameter_names = getattr(result, "parameter_names", None)
    if parameter_names is None:
        return [
            str(column)
            for column in posterior_data.columns
            if str(column) not in {"x", "weights", "loglikelihood"}
        ]
    return [str(name) for name in parameter_names if str(name) in posterior_data.columns]


def _format_corner_label(parameter_name: str, parameter_name_map: dict[str, str] | None) -> str:
    actual_to_base = _invert_parameter_name_map(parameter_name_map)
    base_name = actual_to_base.get(parameter_name, parameter_name)
    return CORNER_LABELS.get(base_name, base_name)


def _apply_corner_teff_formatting(
    figure,
    *,
    parameter_names: list[str],
    labels: list[str],
    samples: np.ndarray,
    weights: np.ndarray | None,
    parameter_name_map: dict[str, str] | None,
    title_kwargs: dict[str, Any] | None = None,
) -> None:
    axes = _reshape_corner_axes(figure, len(parameter_names))
    if axes is None:
        return

    actual_to_base = _invert_parameter_name_map(parameter_name_map)
    for index, parameter_name in enumerate(parameter_names):
        base_name = actual_to_base.get(parameter_name, parameter_name)
        if not _is_teff_parameter(base_name, parameter_name):
            continue

        axes[index, index].set_title(
            _format_teff_corner_title(labels[index], samples[:, index], weights),
            **(title_kwargs or {}),
        )
        for row in range(index, len(parameter_names)):
            _set_teff_axis_formatter(axes[row, index].xaxis)
        for column in range(index + 1):
            _set_teff_axis_formatter(axes[index, column].yaxis)


def _reshape_corner_axes(figure, dimension_count: int) -> np.ndarray | None:
    if dimension_count <= 0:
        return None

    axes = np.asarray(figure.axes, dtype=object)
    expected_count = dimension_count * dimension_count
    if axes.size < expected_count:
        return None
    return axes[:expected_count].reshape(dimension_count, dimension_count)


def _is_teff_parameter(base_name: str, parameter_name: str) -> bool:
    return base_name == "teff" or parameter_name == "teff" or parameter_name.startswith("teff_")


def _format_teff_corner_title(label: str, samples: np.ndarray, weights: np.ndarray | None) -> str:
    q_lo, q_mid, q_hi = _quantile(samples, [0.16, 0.5, 0.84], weights=weights)
    median = int(round(q_mid))
    lower = int(round(q_mid - q_lo))
    upper = int(round(q_hi - q_mid))
    return f"{label} = ${{{median}}}_{{-{lower}}}^{{+{upper}}}$"


def _set_teff_axis_formatter(axis) -> None:
    axis.set_major_formatter(FuncFormatter(_format_teff_tick_value))
    axis.offsetText.set_visible(False)


def _format_teff_tick_value(value: float, _position: int) -> str:
    if not np.isfinite(value):
        return str(value)
    return str(int(round(value)))


def _invert_parameter_name_map(parameter_name_map: dict[str, str] | None) -> dict[str, str]:
    if parameter_name_map is None:
        return {}
    return {actual: base for base, actual in parameter_name_map.items()}


def _quantile(
    samples: np.ndarray,
    q: list[float],
    *,
    weights: np.ndarray | None = None,
) -> np.ndarray:
    values = np.asarray(samples, dtype=float)
    quantiles = np.asarray(q, dtype=float)
    if weights is None:
        return np.percentile(values, 100.0 * quantiles)

    sample_weights = np.asarray(weights, dtype=float)
    if values.shape[0] != sample_weights.shape[0]:
        raise ValueError("Dimension mismatch: len(weights) != len(samples)")
    order = np.argsort(values)
    sorted_values = values[order]
    sorted_weights = sample_weights[order]
    cdf = np.cumsum(sorted_weights)[:-1]
    cdf /= cdf[-1]
    cdf = np.append(0.0, cdf)
    return np.asarray(np.interp(quantiles, cdf, sorted_values), dtype=float)


def _load_corner_module():
    try:
        return importlib.import_module("corner")
    except ImportError:
        return None


def _should_plot_corner_contours(samples: np.ndarray, weights: np.ndarray | None) -> bool:
    sample_count = int(samples.shape[0])
    effective_sample_size = _effective_sample_size(weights, sample_count)
    if sample_count < CORNER_MIN_SAMPLES_FOR_CONTOURS or effective_sample_size < CORNER_MIN_SAMPLES_FOR_CONTOURS:
        LOGGER.info(
            "Disabling corner contours because the posterior sample count is low "
            "(samples=%d, effective_samples=%.1f).",
            sample_count,
            effective_sample_size,
        )
        return False
    return True


def _effective_sample_size(weights: np.ndarray | None, sample_count: int) -> float:
    if weights is None:
        return float(sample_count)

    total_weight = float(np.sum(weights))
    if total_weight <= 0.0:
        return 0.0

    sum_weight_sq = float(np.sum(np.square(weights)))
    if sum_weight_sq <= 0.0:
        return 0.0
    return (total_weight * total_weight) / sum_weight_sq


def _order_summary_parameters(
    summary: dict[str, Any],
    *,
    grid_parameter_order: list[str],
    parameter_name_map: dict[str, str] | None,
) -> dict[str, Any]:
    ordered_keys: list[str] = []
    resolved_parameter_name_map = parameter_name_map or {}

    for grid_parameter in grid_parameter_order:
        summary_key = resolved_parameter_name_map.get(grid_parameter, grid_parameter)
        if summary_key in summary and summary_key not in ordered_keys:
            ordered_keys.append(summary_key)

    for key in summary:
        if key not in ordered_keys:
            ordered_keys.append(key)

    return {key: summary[key] for key in ordered_keys}


def _get_grid_parameter_order(grid) -> list[str]:
    if grid is None:
        return []

    grid_param_names = getattr(grid, "param_names", None)
    if grid_param_names is not None:
        return [str(name) for name in grid_param_names]

    meta_grid = getattr(grid, "meta_grid", None)
    if meta_grid is None:
        return []

    if hasattr(meta_grid, "get"):
        raw_parameters = meta_grid.get("parameters")
        if isinstance(raw_parameters, (list, tuple)):
            return [str(name) for name in raw_parameters]
        if raw_parameters is not None:
            return [str(raw_parameters)]

    return []


def _build_derived_summary(
    *,
    median: dict[str, object],
    mean: dict[str, object],
    parameter_name_map: dict[str, str] | None,
) -> dict[str, object]:
    mapping = parameter_name_map or {}
    logg_key = mapping.get("logg", "logg")
    if logg_key in median or logg_key in mean or "logg" in median or "logg" in mean:
        return {}

    teff_key = mapping.get("teff", "teff")
    mh_key = mapping.get("mh", "mh")
    if teff_key not in median or mh_key not in median or teff_key not in mean or mh_key not in mean:
        return {}

    median_logg = _predict_scalar_logg(float(median[teff_key]), float(median[mh_key]))
    mean_logg = _predict_scalar_logg(float(mean[teff_key]), float(mean[mh_key]))
    median["logg"] = median_logg
    mean["logg"] = mean_logg
    return {
        "logg_source": "yyisochrones_galactic_center",
        "logg_from_grid": False,
    }


def _predict_scalar_logg(teff: float, mh: float) -> float:
    try:
        return float(np.asarray(predict_logg_per_slice(teff, mh, branch="giant")).reshape(-1)[0])
    except ValueError:
        LOGGER.warning(
            "Falling back to global polynomial logg predictor for Teff=%.1f, [M/H]=%.3f",
            teff,
            mh,
        )
        return float(np.asarray(predict_logg(teff, mh)).reshape(-1)[0])
