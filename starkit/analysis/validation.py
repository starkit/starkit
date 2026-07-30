from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
from pathlib import Path
from typing import Any

from .config import (
    AnalysisConfig,
    PriorConfig,
    absolute_user_path,
    format_missing_path_message,
    load_config_data,
)
from starkit.gridkit.inspector import GridInspection, ParameterInspection, inspect_grid

ALLOWED_WAVELENGTH_TYPES = {"air", "vacuum"}
ALLOWED_INPUT_FORMATS = {"ascii_columns", "npy_columns"}
ALLOWED_LIKELIHOODS = {"chi2", "l1"}
ALLOWED_NORMALIZE_PARTS = {"segments", "autodetect"}
OUTPUT_FILENAME_KEYS = {
    "bestfit_spectrum",
    "corner_plot",
    "hdf5",
    "masked_observed",
    "plot",
    "run_summary",
}


@dataclass(frozen=True)
class ValidationIssue:
    level: str
    code: str
    message: str


@dataclass(frozen=True)
class ValidationReport:
    path: Path
    issues: list[ValidationIssue]
    deep_checks_run: bool = False

    @property
    def ok(self) -> bool:
        return not any(issue.level == "error" for issue in self.issues)

    @property
    def errors(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.level == "error"]

    @property
    def warnings(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.level == "warning"]

    @property
    def info(self) -> list[ValidationIssue]:
        return [issue for issue in self.issues if issue.level == "info"]

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "ok": self.ok,
            "deep_checks_run": self.deep_checks_run,
            "errors": [asdict(issue) for issue in self.errors],
            "warnings": [asdict(issue) for issue in self.warnings],
            "info": [asdict(issue) for issue in self.info],
        }


def verify_config(path: str | Path, *, deep: bool = True) -> ValidationReport:
    config_path = absolute_user_path(path)
    issues: list[ValidationIssue] = []
    deep_checks_run = False

    raw = _load_raw_json(config_path, issues)
    if raw is None:
        return ValidationReport(path=config_path, issues=issues, deep_checks_run=deep_checks_run)
    if not isinstance(raw, dict):
        _error(issues, "json.root_type", "Top-level JSON value must be an object.")
        return ValidationReport(path=config_path, issues=issues, deep_checks_run=deep_checks_run)

    config = _load_typed_config(config_path, raw, issues)
    if config is None:
        return ValidationReport(path=config_path, issues=issues, deep_checks_run=deep_checks_run)

    _validate_basic_config(config, issues)
    _validate_prior_shapes(config, issues)
    _validate_component_priors(config, issues)
    _validate_path_references(config, issues)

    inspection: GridInspection | None = None
    if config.run.grid_path.exists():
        inspection = _inspect_grid(config.run.grid_path, issues)
        if inspection is not None:
            _validate_priors_against_grid(config, inspection, issues)

    if deep and not any(issue.level == "error" for issue in issues):
        deep_checks_run = True
        _run_deep_preflight(config, issues)

    if deep_checks_run and not any(issue.level == "error" for issue in issues):
        _info(issues, "deep.ok", "Deep preflight completed successfully.")

    return ValidationReport(path=config_path, issues=issues, deep_checks_run=deep_checks_run)


def format_report(report: ValidationReport) -> str:
    status = "OK" if report.ok else "FAILED"
    lines = [
        f"Status: {status}",
        f"Config: {report.path}",
        f"Checks: {'deep' if report.deep_checks_run else 'static'}",
    ]

    if report.errors:
        lines.append("")
        lines.append("Errors:")
        lines.extend(f"- [{issue.code}] {issue.message}" for issue in report.errors)

    if report.warnings:
        lines.append("")
        lines.append("Warnings:")
        lines.extend(f"- [{issue.code}] {issue.message}" for issue in report.warnings)

    if report.info:
        lines.append("")
        lines.append("Info:")
        lines.extend(f"- [{issue.code}] {issue.message}" for issue in report.info)

    return "\n".join(lines)


def _load_raw_json(config_path: Path, issues: list[ValidationIssue]) -> dict[str, Any] | list[Any] | None:
    if not config_path.exists():
        _error(issues, "config.file_missing", f"Config file does not exist: {config_path}")
        return None

    try:
        text = config_path.read_text(encoding="ascii")
    except UnicodeDecodeError as exc:
        _error(issues, "json.encoding", f"Config file must be ASCII-decodable: {exc}")
        return None
    except OSError as exc:
        _error(issues, "config.read", f"Could not read config file {config_path}: {exc}")
        return None

    duplicate_keys: list[str] = []

    def _hook(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        data: dict[str, Any] = {}
        seen: set[str] = set()
        for key, value in pairs:
            if key in seen:
                duplicate_keys.append(key)
            seen.add(key)
            data[key] = value
        return data

    try:
        raw = json.loads(text, object_pairs_hook=_hook)
    except json.JSONDecodeError as exc:
        _error(
            issues,
            "json.syntax",
            f"Invalid JSON syntax at line {exc.lineno}, column {exc.colno}: {exc.msg}",
        )
        return None

    for key in sorted(set(duplicate_keys)):
        _error(
            issues,
            "json.duplicate_key",
            f"Duplicate JSON key {key!r} detected; later values overwrite earlier ones.",
        )

    return raw


def _load_typed_config(
    config_path: Path,
    raw: dict[str, Any],
    issues: list[ValidationIssue],
) -> AnalysisConfig | None:
    try:
        return load_config_data(raw, config_path=config_path)
    except Exception as exc:
        _error(issues, "config.schema", str(exc))
        return None


def _validate_basic_config(config: AnalysisConfig, issues: list[ValidationIssue]) -> None:
    _check_enum(issues, "run.grid_wavelength_type", config.run.grid_wavelength_type, ALLOWED_WAVELENGTH_TYPES)
    _check_enum(issues, "input_spectrum.format", config.input_spectrum.format, ALLOWED_INPUT_FORMATS, required=True)
    _check_enum(issues, "input_spectrum.wavelength_type", config.input_spectrum.wavelength_type, ALLOWED_WAVELENGTH_TYPES)
    _check_enum(issues, "segments.wavelength_type", config.segments.wavelength_type, ALLOWED_WAVELENGTH_TYPES)
    _check_enum(issues, "masking.wavelength_type", config.masking.wavelength_type, ALLOWED_WAVELENGTH_TYPES)
    _check_enum(issues, "outputs.wavelength_type", config.outputs.wavelength_type, ALLOWED_WAVELENGTH_TYPES)
    _check_enum(issues, "model.likelihood", config.model.likelihood, ALLOWED_LIKELIHOODS, required=True)
    _check_enum(issues, "model.normalize_parts", config.model.normalize_parts, ALLOWED_NORMALIZE_PARTS)

    if not _is_integer_like(config.run.mpi_processes):
        _error(issues, "run.mpi_processes", "run.mpi_processes must be an integer.")
    elif int(float(config.run.mpi_processes)) < 1:
        _error(issues, "run.mpi_processes", "run.mpi_processes must be at least 1.")
    elif not config.run.mpi and int(float(config.run.mpi_processes)) != 1:
        _warning(
            issues,
            "run.mpi_processes.ignored",
            f"run.mpi is false, so mpi_processes={config.run.mpi_processes} will be ignored.",
        )

    if config.input_spectrum.format in ALLOWED_INPUT_FORMATS:
        for name in ("wavelength", "flux"):
            if config.input_spectrum.columns.get(name) is None:
                _error(
                    issues,
                    f"input.columns.{name}",
                    f"input_spectrum.columns.{name} must be configured for format {config.input_spectrum.format!r}.",
                )
        for name in ("wavelength", "flux", "uncertainty"):
            value = config.input_spectrum.columns.get(name)
            if value is not None and not _is_integer_like(value):
                _error(
                    issues,
                    f"input.columns.{name}",
                    f"input_spectrum.columns.{name} must be an integer column index.",
                )

    if not _is_integer_like(config.input_spectrum.skiprows):
        _error(issues, "input.skiprows", "input_spectrum.skiprows must be an integer.")
    elif int(float(config.input_spectrum.skiprows)) < 0:
        _error(issues, "input.skiprows", "input_spectrum.skiprows must be zero or positive.")

    if (
        config.input_spectrum.columns.get("uncertainty") is not None
        and config.input_spectrum.uncertainty_path is not None
    ):
        _warning(
            issues,
            "input.uncertainty_path.ignored",
            "input_spectrum.uncertainty_path is ignored because columns.uncertainty is already configured.",
        )
    if config.input_spectrum.uncertainty_column is not None and config.input_spectrum.uncertainty_path is None:
        _warning(
            issues,
            "input.uncertainty_column.ignored",
            "input_spectrum.uncertainty_column is ignored unless uncertainty_path is set.",
        )
    if config.input_spectrum.uncertainty_column is not None and not _is_integer_like(config.input_spectrum.uncertainty_column):
        _error(
            issues,
            "input.uncertainty_column",
            "input_spectrum.uncertainty_column must be an integer column index.",
        )

    has_data_uncertainty = (
        config.input_spectrum.columns.get("uncertainty") is not None
        or config.input_spectrum.uncertainty_path is not None
    )
    if has_data_uncertainty and config.input_spectrum.snr is not None:
        _warning(
            issues,
            "input.snr.ignored",
            "input_spectrum.snr is ignored because an uncertainty source is already configured.",
        )
    if not has_data_uncertainty and config.input_spectrum.snr is None:
        _warning(
            issues,
            "input.uncertainty.missing",
            "No uncertainty column, uncertainty_path, or snr is configured; the observed spectrum will have no uncertainties.",
        )
    if config.input_spectrum.snr is not None:
        if not _is_finite_number(config.input_spectrum.snr):
            _error(issues, "input.snr", "input_spectrum.snr must be numeric.")
        elif float(config.input_spectrum.snr) <= 0:
            _error(issues, "input.snr", "input_spectrum.snr must be positive.")

    if config.segments.mode != "file":
        _warning(
            issues,
            "segments.mode.unsupported",
            f"segments.mode={config.segments.mode!r} is not implemented; the analyser currently behaves as file mode only.",
        )
    if config.segments.path is None and config.segments.resample_to_grid:
        _warning(
            issues,
            "segments.resample.ignored",
            "segments.resample_to_grid is true, but no segments.path is configured, so the full observed wavelength grid will be used.",
        )
    if config.segments.wavelength_reference_path is not None and not config.segments.resample_to_grid:
        _warning(
            issues,
            "segments.reference.ignored",
            "segments.wavelength_reference_path is only used when segments.resample_to_grid is true.",
        )

    normalize_npol = config.model.normalize_npol
    if normalize_npol is not None:
        if not _is_integer_like(normalize_npol):
            _error(issues, "model.normalize_npol", "model.normalize_npol must be an integer.")
        elif int(normalize_npol) < 0:
            _error(issues, "model.normalize_npol", "model.normalize_npol must be zero or positive.")

    if config.model.normalize_parts is not None:
        if normalize_npol is None:
            _error(
                issues,
                "model.normalize_parts.requires_npol",
                "model.normalize_parts requires model.normalize_npol to be set.",
            )
        if config.model.normalize_parts == "segments" and config.segments.path is None:
            _error(
                issues,
                "model.normalize_parts.requires_segments",
                "model.normalize_parts='segments' requires segments.path to be configured.",
            )

    if config.model.instrument_resolution is not None:
        if not _is_finite_number(config.model.instrument_resolution):
            _error(issues, "model.instrument_resolution", "model.instrument_resolution must be numeric.")
        elif float(config.model.instrument_resolution) <= 0:
            _error(issues, "model.instrument_resolution", "model.instrument_resolution must be positive.")

    if not _is_finite_number(config.model.rotational_velocity):
        _error(issues, "model.rotational_velocity", "model.rotational_velocity must be numeric.")
    elif float(config.model.rotational_velocity) < 0:
        _error(issues, "model.rotational_velocity", "model.rotational_velocity must be zero or positive.")
    elif not config.model.rotational_broadening and float(config.model.rotational_velocity) != 0:
        _warning(
            issues,
            "model.rotational_velocity.ignored",
            "model.rotational_velocity is ignored because model.rotational_broadening is false.",
        )

    if not _is_finite_number(config.model.radial_velocity):
        _error(issues, "model.radial_velocity", "model.radial_velocity must be numeric.")
    elif not config.model.doppler_shift and float(config.model.radial_velocity) != 0:
        _warning(
            issues,
            "model.radial_velocity.ignored",
            "model.radial_velocity is ignored because model.doppler_shift is false.",
        )

    unknown_output_keys = sorted(set(config.outputs.filenames) - OUTPUT_FILENAME_KEYS)
    for key in unknown_output_keys:
        _warning(
            issues,
            "outputs.filename_key.unknown",
            f"outputs.filenames[{key!r}] is not a recognized output override key.",
        )

    if not any(
        (
            config.outputs.write_hdf5,
            config.outputs.write_plot,
            config.outputs.write_corner_plot,
            config.outputs.write_bestfit_spectrum,
            config.outputs.write_masked_observed,
            config.outputs.write_run_summary,
        )
    ):
        _warning(issues, "outputs.disabled", "All output writers are disabled; the run would produce no saved outputs.")


def _validate_prior_shapes(config: AnalysisConfig, issues: list[ValidationIssue]) -> None:
    for name, prior in config.priors.items():
        kind = prior.type.lower()
        if kind == "fixed":
            if prior.value is None:
                _error(issues, f"prior.{name}.value", f"Prior {name!r} with type 'fixed' requires 'value'.")
                continue
            if not _is_finite_number(prior.value):
                _error(issues, f"prior.{name}.value", f"Prior {name!r} has a non-finite fixed value.")
            _warn_ignored_prior_fields(issues, name, kind, prior, allowed_fields={"value"})
            continue

        if kind == "uniform":
            if prior.min is None or prior.max is None:
                _error(issues, f"prior.{name}.bounds", f"Prior {name!r} with type 'uniform' requires 'min' and 'max'.")
                continue
            if not _is_finite_number(prior.min) or not _is_finite_number(prior.max):
                _error(issues, f"prior.{name}.bounds", f"Prior {name!r} has non-finite uniform bounds.")
                continue
            if float(prior.min) > float(prior.max):
                _error(
                    issues,
                    f"prior.{name}.bounds",
                    f"Prior {name!r} has min={prior.min} larger than max={prior.max}.",
                )
            elif float(prior.min) == float(prior.max):
                _warning(
                    issues,
                    f"prior.{name}.degenerate_uniform",
                    f"Prior {name!r} has identical uniform bounds; consider a fixed prior instead.",
                )
            _warn_ignored_prior_fields(issues, name, kind, prior, allowed_fields={"min", "max"})
            continue

        if kind == "gaussian":
            if prior.mean is None or prior.sigma is None:
                _error(
                    issues,
                    f"prior.{name}.gaussian",
                    f"Prior {name!r} with type 'gaussian' requires 'mean' and 'sigma'.",
                )
                continue
            if not _is_finite_number(prior.mean) or not _is_finite_number(prior.sigma):
                _error(issues, f"prior.{name}.gaussian", f"Prior {name!r} has non-finite Gaussian parameters.")
                continue
            if float(prior.sigma) <= 0:
                _error(issues, f"prior.{name}.gaussian", f"Prior {name!r} must have sigma > 0.")
            _warn_ignored_prior_fields(issues, name, kind, prior, allowed_fields={"mean", "sigma"})
            continue

        _error(
            issues,
            f"prior.{name}.type",
            f"Prior {name!r} has unsupported type {prior.type!r}; expected fixed, uniform, or gaussian.",
        )


def _validate_component_priors(config: AnalysisConfig, issues: list[ValidationIssue]) -> None:
    if "vrot" in config.priors and not config.model.rotational_broadening:
        _error(
            issues,
            "model.vrot.disabled",
            "Prior 'vrot' is configured, but model.rotational_broadening is false.",
        )
    if config.model.rotational_broadening and "vrot" not in config.priors:
        _warning(
            issues,
            "model.vrot.fixed",
            f"model.rotational_broadening is true, but no 'vrot' prior is configured, so vrot stays fixed at {config.model.rotational_velocity}.",
        )

    if "vrad" in config.priors and not config.model.doppler_shift:
        _error(
            issues,
            "model.vrad.disabled",
            "Prior 'vrad' is configured, but model.doppler_shift is false.",
        )
    if config.model.doppler_shift and "vrad" not in config.priors:
        _warning(
            issues,
            "model.vrad.fixed",
            f"model.doppler_shift is true, but no 'vrad' prior is configured, so vrad stays fixed at {config.model.radial_velocity}.",
        )

    if "R" in config.priors and config.model.instrument_resolution is None:
        _error(
            issues,
            "model.R.disabled",
            "Prior 'R' is configured, but model.instrument_resolution is null.",
        )
    if config.model.instrument_resolution is not None and "R" not in config.priors:
        _warning(
            issues,
            "model.R.fixed",
            f"model.instrument_resolution is set, but no 'R' prior is configured, so R stays fixed at {config.model.instrument_resolution}.",
        )

    if "vrot" in config.priors:
        _check_non_negative_prior("vrot", config.priors["vrot"], issues)
    if "R" in config.priors:
        _check_strictly_positive_prior("R", config.priors["R"], issues)


def _validate_path_references(config: AnalysisConfig, issues: list[ValidationIssue]) -> None:
    _require_existing_path(issues, "run.grid_path", config.run.grid_path, config.source_path)
    _require_existing_path(issues, "input_spectrum.path", config.input_spectrum.path, config.source_path)

    optional_paths = [
        ("input_spectrum.uncertainty_path", config.input_spectrum.uncertainty_path),
        ("segments.path", config.segments.path),
        ("segments.wavelength_reference_path", config.segments.wavelength_reference_path),
        ("masking.include_regions_path", config.masking.include_regions_path),
        ("masking.exclude_regions_path", config.masking.exclude_regions_path),
        ("masking.regions_json_path", config.masking.regions_json_path),
    ]
    for field_name, path in optional_paths:
        if path is not None:
            _require_existing_path(issues, field_name, path, config.source_path)


def _inspect_grid(grid_path: Path, issues: list[ValidationIssue]) -> GridInspection | None:
    try:
        inspection = inspect_grid(grid_path)
    except Exception as exc:
        _error(issues, "grid.inspect", f"Could not inspect grid {grid_path}: {exc}")
        return None

    for warning in inspection.warnings:
        _warning(issues, "grid.warning", warning)
    return inspection


def _validate_priors_against_grid(
    config: AnalysisConfig,
    inspection: GridInspection,
    issues: list[ValidationIssue],
) -> None:
    parameter_lookup = {parameter.name: parameter for parameter in inspection.parameters}
    for name, prior in config.priors.items():
        if name in parameter_lookup:
            _check_grid_parameter_prior(name, prior, parameter_lookup[name], issues)

    resolution = inspection.resolution
    if resolution is None:
        return

    if (
        config.model.instrument_resolution is not None
        and _is_finite_number(config.model.instrument_resolution)
        and float(config.model.instrument_resolution) >= float(resolution)
    ):
        _error(
            issues,
            "model.instrument_resolution.grid_limit",
            "model.instrument_resolution must be smaller than the grid resolution "
            f"(requested {config.model.instrument_resolution}, grid {resolution}).",
        )

    resolution_prior = config.priors.get("R")
    if resolution_prior is None:
        return

    kind = resolution_prior.type.lower()
    if (
        kind == "uniform"
        and resolution_prior.max is not None
        and _is_finite_number(resolution_prior.max)
        and float(resolution_prior.max) >= float(resolution)
    ):
        _error(
            issues,
            "prior.R.above_grid_resolution",
            "Uniform prior for R must have max smaller than the grid resolution "
            f"(prior max {resolution_prior.max}, grid {resolution}).",
        )
    if (
        kind == "fixed"
        and resolution_prior.value is not None
        and _is_finite_number(resolution_prior.value)
        and float(resolution_prior.value) >= float(resolution)
    ):
        _error(
            issues,
            "prior.R.above_grid_resolution",
            "Fixed prior for R must be smaller than the grid resolution "
            f"(prior value {resolution_prior.value}, grid {resolution}).",
        )
    if (
        kind == "gaussian"
        and resolution_prior.mean is not None
        and resolution_prior.sigma is not None
        and _is_finite_number(resolution_prior.mean)
        and _is_finite_number(resolution_prior.sigma)
    ):
        upper_3sigma = float(resolution_prior.mean) + 3.0 * float(resolution_prior.sigma)
        if upper_3sigma >= float(resolution):
            _warning(
                issues,
                "prior.R.gaussian_grid_edge",
                "Gaussian prior for R extends close to or beyond the grid resolution "
                f"(mean + 3 sigma = {upper_3sigma}, grid {resolution}).",
            )


def _check_grid_parameter_prior(
    name: str,
    prior: PriorConfig,
    parameter: ParameterInspection,
    issues: list[ValidationIssue],
) -> None:
    lower = float(parameter.minimum)
    upper = float(parameter.maximum)
    kind = prior.type.lower()

    if kind == "fixed" and prior.value is not None and _is_finite_number(prior.value):
        value = float(prior.value)
        if value < lower or value > upper:
            _error(
                issues,
                f"prior.{name}.outside_grid",
                f"Fixed prior for {name!r} is outside the grid range [{lower}, {upper}] (value {value}).",
            )
        return

    if (
        kind == "uniform"
        and prior.min is not None
        and prior.max is not None
        and _is_finite_number(prior.min)
        and _is_finite_number(prior.max)
    ):
        if float(prior.min) < lower or float(prior.max) > upper:
            _error(
                issues,
                f"prior.{name}.outside_grid",
                f"Uniform prior for {name!r} extends outside the grid range [{lower}, {upper}] "
                f"(requested [{prior.min}, {prior.max}]).",
            )
        return

    if (
        kind == "gaussian"
        and prior.mean is not None
        and prior.sigma is not None
        and _is_finite_number(prior.mean)
        and _is_finite_number(prior.sigma)
    ):
        mean = float(prior.mean)
        sigma = float(prior.sigma)
        if mean < lower or mean > upper:
            _warning(
                issues,
                f"prior.{name}.gaussian_mean_outside_grid",
                f"Gaussian prior for {name!r} has mean {mean} outside the grid range [{lower}, {upper}].",
            )
            return
        lower_3sigma = mean - 3.0 * sigma
        upper_3sigma = mean + 3.0 * sigma
        if lower_3sigma < lower or upper_3sigma > upper:
            _warning(
                issues,
                f"prior.{name}.gaussian_grid_edge",
                f"Gaussian prior for {name!r} extends close to or beyond the grid range "
                f"([{lower_3sigma}, {upper_3sigma}] vs. [{lower}, {upper}]).",
            )


def _run_deep_preflight(config: AnalysisConfig, issues: list[ValidationIssue]) -> None:
    from .pipeline import build_model, get_grid_sampling, load_analysis_grid
    from .spectrum import (
        apply_masking,
        autodetect_normalization_segment_bounds,
        build_segment_normalization_parts,
        build_observed_spectrum,
        load_input_arrays,
        prepare_analysis_arrays,
    )

    try:
        raw_arrays = load_input_arrays(config.input_spectrum)
    except Exception as exc:
        _error(issues, "deep.input", f"Could not load input spectrum: {exc}")
        return

    try:
        grid = load_analysis_grid(config)
    except Exception as exc:
        _error(issues, "deep.grid", f"Could not load analysis grid: {exc}")
        return

    grid_wavelength_type = getattr(grid, "wavelength_type", None)
    grid_resolution, grid_sampling = get_grid_sampling(grid)

    try:
        analysis_arrays = prepare_analysis_arrays(
            raw_arrays=raw_arrays,
            input_config=config.input_spectrum,
            segments_config=config.segments,
            grid_resolution=grid_resolution,
            grid_sampling=grid_sampling,
            grid_reference_wavelength=grid.wavelength.value,
            grid_wavelength_type=grid_wavelength_type,
        )
    except Exception as exc:
        _error(issues, "deep.segments", f"Could not prepare analysis wavelength arrays: {exc}")
        return

    analysis_wavelength_type = grid_wavelength_type or config.input_spectrum.wavelength_type
    try:
        masked_arrays = apply_masking(
            analysis_arrays.wavelength,
            analysis_arrays.flux,
            analysis_arrays.uncertainty,
            config.masking,
            wavelength_unit=config.input_spectrum.wavelength_unit,
            default_wavelength_type=config.input_spectrum.wavelength_type,
            target_wavelength_type=analysis_wavelength_type,
            segment_bounds=analysis_arrays.segment_bounds,
        )
    except Exception as exc:
        _error(issues, "deep.masking", f"Could not apply masking: {exc}")
        return

    try:
        observed = build_observed_spectrum(masked_arrays, config.input_spectrum)
    except Exception as exc:
        _error(issues, "deep.observed", f"Could not build the observed spectrum object: {exc}")
        return

    normalization_parts = None
    if config.model.normalize_parts in {"segments", "autodetect"} and config.model.normalize_npol is not None:
        try:
            segment_bounds = masked_arrays.segment_bounds
            if config.model.normalize_parts == "autodetect":
                segment_bounds = autodetect_normalization_segment_bounds(analysis_arrays.wavelength)
            normalization_parts = build_segment_normalization_parts(
                observed,
                segment_bounds,
                config.model.normalize_npol,
            )
            if config.model.normalize_parts == "autodetect" and len(normalization_parts) == 1:
                _warning(
                    issues,
                    "deep.normalize_parts.single_part",
                    "model.normalize_parts='autodetect' found no large gaps, so normalization will behave like one global part.",
                )
        except Exception as exc:
            _error(issues, "deep.normalize_parts", f"Could not build normalization parts: {exc}")
            return

    try:
        build_model(config, grid, observed, normalization_parts=normalization_parts)
    except Exception as exc:
        _error(issues, "deep.model", f"Could not build the fit model: {exc}")


def _warn_ignored_prior_fields(
    issues: list[ValidationIssue],
    name: str,
    kind: str,
    prior: PriorConfig,
    *,
    allowed_fields: set[str],
) -> None:
    field_names = {
        "value": prior.value,
        "min": prior.min,
        "max": prior.max,
        "mean": prior.mean,
        "sigma": prior.sigma,
    }
    ignored = sorted(field_name for field_name, value in field_names.items() if value is not None and field_name not in allowed_fields)
    if ignored:
        _warning(
            issues,
            f"prior.{name}.ignored_fields",
            f"Prior {name!r} with type {kind!r} also sets ignored fields: {', '.join(ignored)}.",
        )


def _check_non_negative_prior(name: str, prior: PriorConfig, issues: list[ValidationIssue]) -> None:
    kind = prior.type.lower()
    if kind == "fixed" and prior.value is not None and _is_finite_number(prior.value) and float(prior.value) < 0:
        _error(issues, f"prior.{name}.negative", f"Prior {name!r} must be zero or positive.")
    if kind == "uniform" and prior.min is not None and _is_finite_number(prior.min) and float(prior.min) < 0:
        _error(issues, f"prior.{name}.negative", f"Prior {name!r} must have min >= 0.")
    if (
        kind == "gaussian"
        and prior.mean is not None
        and prior.sigma is not None
        and _is_finite_number(prior.mean)
        and _is_finite_number(prior.sigma)
    ):
        lower_3sigma = float(prior.mean) - 3.0 * float(prior.sigma)
        if lower_3sigma < 0:
            _warning(
                issues,
                f"prior.{name}.negative_tail",
                f"Gaussian prior {name!r} extends below zero (mean - 3 sigma = {lower_3sigma}).",
            )


def _check_strictly_positive_prior(name: str, prior: PriorConfig, issues: list[ValidationIssue]) -> None:
    kind = prior.type.lower()
    if kind == "fixed" and prior.value is not None and _is_finite_number(prior.value) and float(prior.value) <= 0:
        _error(issues, f"prior.{name}.non_positive", f"Prior {name!r} must be > 0.")
    if kind == "uniform" and prior.min is not None and _is_finite_number(prior.min) and float(prior.min) <= 0:
        _error(issues, f"prior.{name}.non_positive", f"Prior {name!r} must have min > 0.")
    if (
        kind == "gaussian"
        and prior.mean is not None
        and prior.sigma is not None
        and _is_finite_number(prior.mean)
        and _is_finite_number(prior.sigma)
    ):
        lower_3sigma = float(prior.mean) - 3.0 * float(prior.sigma)
        if lower_3sigma <= 0:
            _warning(
                issues,
                f"prior.{name}.non_positive_tail",
                f"Gaussian prior {name!r} extends to non-positive values (mean - 3 sigma = {lower_3sigma}).",
            )


def _require_existing_path(
    issues: list[ValidationIssue],
    field_name: str,
    path: Path,
    config_path: Path | None,
) -> None:
    if not path.exists():
        _error(
            issues,
            f"{field_name}.missing",
            format_missing_path_message(field_name, path, config_path=config_path),
        )


def _check_enum(
    issues: list[ValidationIssue],
    field_name: str,
    value: str | None,
    allowed: set[str],
    *,
    required: bool = False,
) -> None:
    if value is None:
        if required:
            _error(issues, field_name, f"{field_name} is required.")
        return
    if value not in allowed:
        allowed_text = ", ".join(sorted(allowed))
        _error(issues, field_name, f"{field_name} must be one of: {allowed_text}. Got {value!r}.")


def _is_finite_number(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def _is_integer_like(value: Any) -> bool:
    if isinstance(value, bool):
        return False
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return False
    return math.isfinite(numeric) and numeric.is_integer()


def _error(issues: list[ValidationIssue], code: str, message: str) -> None:
    issues.append(ValidationIssue(level="error", code=code, message=message))


def _warning(issues: list[ValidationIssue], code: str, message: str) -> None:
    issues.append(ValidationIssue(level="warning", code=code, message=message))


def _info(issues: list[ValidationIssue], code: str, message: str) -> None:
    issues.append(ValidationIssue(level="info", code=code, message=message))
