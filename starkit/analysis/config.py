from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class PriorConfig:
    type: str
    value: float | None = None
    min: float | None = None
    max: float | None = None
    mean: float | None = None
    sigma: float | None = None


@dataclass(frozen=True)
class RunConfig:
    name: str
    grid_path: Path
    grid_wavelength_type: str | None = None
    mpi: bool = True
    mpi_processes: int = 1
    work_dir: Path | None = None
    clean_up: bool = True


@dataclass(frozen=True)
class InputSpectrumConfig:
    path: Path
    format: str
    wavelength_unit: str
    flux_unit: str
    wavelength_type: str | None = None
    columns: dict[str, int | None] = field(default_factory=dict)
    skiprows: int = 0
    uncertainty_path: Path | None = None
    uncertainty_column: int | None = None
    snr: float | None = None


@dataclass(frozen=True)
class SegmentsConfig:
    path: Path | None = None
    wavelength_reference_path: Path | None = None
    wavelength_type: str | None = None
    mode: str = "file"
    resample_to_grid: bool = False


@dataclass(frozen=True)
class MaskingConfig:
    wavelength_type: str | None = None
    include_regions_path: Path | None = None
    exclude_regions_path: Path | None = None
    regions_json_path: Path | None = None


@dataclass(frozen=True)
class ModelConfig:
    likelihood: str = "chi2"
    crop_to_observed_range: bool = True
    interpolate: bool = True
    normalize_npol: int | None = None
    normalize_parts: str | None = None
    instrument_resolution: float | None = None
    rotational_broadening: bool = False
    rotational_velocity: float = 0.0
    doppler_shift: bool = False
    radial_velocity: float = 0.0


@dataclass(frozen=True)
class OutputConfig:
    directory: Path
    wavelength_type: str | None = None
    write_hdf5: bool = True
    write_plot: bool = True
    write_corner_plot: bool = True
    write_bestfit_spectrum: bool = True
    write_masked_observed: bool = True
    write_run_summary: bool = True
    filenames: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class AnalysisConfig:
    run: RunConfig
    input_spectrum: InputSpectrumConfig
    segments: SegmentsConfig
    masking: MaskingConfig
    model: ModelConfig
    priors: dict[str, PriorConfig]
    outputs: OutputConfig
    source_path: Path | None = None


ROOT_KEYS = {
    "run",
    "input_spectrum",
    "segments",
    "masking",
    "model",
    "priors",
    "outputs",
}
RUN_KEYS = {
    "name",
    "grid_path",
    "grid_wavelength_type",
    "mpi",
    "mpi_processes",
    "work_dir",
    "clean_up",
}
INPUT_SPECTRUM_KEYS = {
    "path",
    "format",
    "wavelength_unit",
    "flux_unit",
    "wavelength_type",
    "columns",
    "skiprows",
    "uncertainty_path",
    "uncertainty_column",
    "snr",
}
INPUT_COLUMN_KEYS = {"wavelength", "flux", "uncertainty"}
SEGMENTS_KEYS = {
    "path",
    "wavelength_reference_path",
    "wavelength_type",
    "mode",
    "resample_to_grid",
}
MASKING_KEYS = {
    "wavelength_type",
    "include_regions_path",
    "exclude_regions_path",
    "regions_json_path",
}
MODEL_KEYS = {
    "likelihood",
    "crop_to_observed_range",
    "interpolate",
    "normalize_npol",
    "normalize_parts",
    "instrument_resolution",
    "rotational_broadening",
    "rotational_velocity",
    "doppler_shift",
    "radial_velocity",
}
PRIOR_KEYS = {"type", "value", "min", "max", "mean", "sigma"}
OUTPUT_KEYS = {
    "directory",
    "wavelength_type",
    "write_hdf5",
    "write_plot",
    "write_corner_plot",
    "write_bestfit_spectrum",
    "write_masked_observed",
    "write_run_summary",
    "filenames",
}


def _as_path(value: str | None, *, config_dir: Path) -> Path | None:
    if value is None:
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = config_dir / path
    return path


def absolute_user_path(path: str | Path) -> Path:
    candidate = Path(path).expanduser()
    if candidate.is_absolute():
        return candidate
    return Path.cwd() / candidate


def _require_keys(raw: dict[str, Any], section: str, keys: tuple[str, ...]) -> None:
    missing = [key for key in keys if key not in raw]
    if missing:
        joined = ", ".join(missing)
        raise ValueError(f"Missing required keys in '{section}': {joined}")


def _require_mapping(value: Any, section: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"Section '{section}' must be a JSON object.")
    return value


def _require_known_keys(raw: dict[str, Any], section: str, allowed: set[str]) -> None:
    unknown = sorted(set(raw) - allowed)
    if unknown:
        raise ValueError(f"Unknown keys in '{section}': {', '.join(unknown)}")


def _require_string(raw: dict[str, Any], section: str, key: str) -> str:
    value = raw[key]
    if not isinstance(value, str):
        raise ValueError(f"{section}.{key} must be a string.")
    return value


def _optional_string(raw: dict[str, Any], section: str, key: str) -> str | None:
    value = raw.get(key)
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"{section}.{key} must be a string or null.")
    return value


def _require_bool(raw: dict[str, Any], section: str, key: str, default: bool) -> bool:
    value = raw.get(key, default)
    if not isinstance(value, bool):
        raise ValueError(f"{section}.{key} must be a boolean.")
    return value


def _is_number(value: Any) -> bool:
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(float(value))


def _require_integer_like(raw: dict[str, Any], section: str, key: str, default: int | None = None) -> int:
    if key not in raw:
        if default is None:
            raise ValueError(f"{section}.{key} is required.")
        return default

    value = raw[key]
    if not _is_number(value) or not float(value).is_integer():
        raise ValueError(f"{section}.{key} must be an integer.")
    return int(value)


def _optional_integer_like(raw: dict[str, Any], section: str, key: str) -> int | None:
    if key not in raw or raw[key] is None:
        return None
    value = raw[key]
    if not _is_number(value) or not float(value).is_integer():
        raise ValueError(f"{section}.{key} must be an integer or null.")
    return int(value)


def _optional_number(raw: dict[str, Any], section: str, key: str) -> float | None:
    if key not in raw or raw[key] is None:
        return None
    value = raw[key]
    if not _is_number(value):
        raise ValueError(f"{section}.{key} must be numeric or null.")
    return float(value)


def _require_path(raw: dict[str, Any], section: str, key: str, *, config_dir: Path) -> Path:
    value = _require_string(raw, section, key)
    return _as_path(value, config_dir=config_dir) or Path(".")


def _optional_path(raw: dict[str, Any], section: str, key: str, *, config_dir: Path) -> Path | None:
    value = _optional_string(raw, section, key)
    return _as_path(value, config_dir=config_dir)


def _load_json_with_duplicate_detection(config_path: Path) -> Any:
    text = config_path.read_text(encoding="ascii")
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
        raise ValueError(
            f"Invalid JSON syntax at line {exc.lineno}, column {exc.colno}: {exc.msg}"
        ) from exc

    if duplicate_keys:
        duplicate_text = ", ".join(sorted(set(duplicate_keys)))
        raise ValueError(f"Duplicate JSON key(s) detected: {duplicate_text}")

    return raw


def format_path_resolution_hint(config_path: Path | None) -> str:
    if config_path is None:
        return "Relative paths in the config are resolved from the config file directory, not the current working directory."
    return (
        f"Config file: {config_path}. "
        f"Relative paths in the config are resolved from {config_path.parent}, not the current working directory."
    )


def format_missing_path_message(field_name: str, path: Path, *, config_path: Path | None) -> str:
    return (
        f"{field_name} does not exist: {path}. "
        f"{format_path_resolution_hint(config_path)}"
    )


def load_config_data(raw: Any, *, config_path: Path) -> AnalysisConfig:
    config_dir = config_path.parent
    root = _require_mapping(raw, "root")
    _require_known_keys(root, "root", ROOT_KEYS)
    _require_keys(root, "root", ("run", "input_spectrum", "model", "priors", "outputs"))

    run_raw = _require_mapping(root["run"], "run")
    _require_known_keys(run_raw, "run", RUN_KEYS)
    _require_keys(run_raw, "run", ("name", "grid_path"))
    run = RunConfig(
        name=_require_string(run_raw, "run", "name"),
        grid_path=_require_path(run_raw, "run", "grid_path", config_dir=config_dir),
        grid_wavelength_type=_optional_string(run_raw, "run", "grid_wavelength_type"),
        mpi=_require_bool(run_raw, "run", "mpi", True),
        mpi_processes=_require_integer_like(run_raw, "run", "mpi_processes", 1),
        work_dir=_optional_path(run_raw, "run", "work_dir", config_dir=config_dir),
        clean_up=_require_bool(run_raw, "run", "clean_up", True),
    )
    if run.mpi_processes < 1:
        raise ValueError("run.mpi_processes must be at least 1")

    input_raw = _require_mapping(root["input_spectrum"], "input_spectrum")
    _require_known_keys(input_raw, "input_spectrum", INPUT_SPECTRUM_KEYS)
    _require_keys(
        input_raw,
        "input_spectrum",
        ("path", "format", "wavelength_unit", "flux_unit"),
    )
    columns_raw = input_raw.get("columns", {})
    columns = _require_mapping(columns_raw, "input_spectrum.columns")
    _require_known_keys(columns, "input_spectrum.columns", INPUT_COLUMN_KEYS)
    input_spectrum = InputSpectrumConfig(
        path=_require_path(input_raw, "input_spectrum", "path", config_dir=config_dir),
        format=_require_string(input_raw, "input_spectrum", "format"),
        wavelength_unit=_require_string(input_raw, "input_spectrum", "wavelength_unit"),
        flux_unit=_require_string(input_raw, "input_spectrum", "flux_unit"),
        wavelength_type=_optional_string(input_raw, "input_spectrum", "wavelength_type"),
        columns={
            name: _optional_integer_like(columns, "input_spectrum.columns", name)
            for name in INPUT_COLUMN_KEYS
            if name in columns
        },
        skiprows=_require_integer_like(input_raw, "input_spectrum", "skiprows", 0),
        uncertainty_path=_optional_path(input_raw, "input_spectrum", "uncertainty_path", config_dir=config_dir),
        uncertainty_column=_optional_integer_like(input_raw, "input_spectrum", "uncertainty_column"),
        snr=_optional_number(input_raw, "input_spectrum", "snr"),
    )

    segments_raw = _require_mapping(root.get("segments", {}), "segments")
    _require_known_keys(segments_raw, "segments", SEGMENTS_KEYS)
    segments = SegmentsConfig(
        path=_optional_path(segments_raw, "segments", "path", config_dir=config_dir),
        wavelength_reference_path=_optional_path(
            segments_raw,
            "segments",
            "wavelength_reference_path",
            config_dir=config_dir,
        ),
        wavelength_type=_optional_string(segments_raw, "segments", "wavelength_type"),
        mode=_optional_string(segments_raw, "segments", "mode") or "file",
        resample_to_grid=_require_bool(segments_raw, "segments", "resample_to_grid", False),
    )

    masking_raw = _require_mapping(root.get("masking", {}), "masking")
    _require_known_keys(masking_raw, "masking", MASKING_KEYS)
    masking = MaskingConfig(
        wavelength_type=_optional_string(masking_raw, "masking", "wavelength_type"),
        include_regions_path=_optional_path(masking_raw, "masking", "include_regions_path", config_dir=config_dir),
        exclude_regions_path=_optional_path(masking_raw, "masking", "exclude_regions_path", config_dir=config_dir),
        regions_json_path=_optional_path(masking_raw, "masking", "regions_json_path", config_dir=config_dir),
    )

    model_raw = _require_mapping(root["model"], "model")
    _require_known_keys(model_raw, "model", MODEL_KEYS)
    model = ModelConfig(
        likelihood=_optional_string(model_raw, "model", "likelihood") or "chi2",
        crop_to_observed_range=_require_bool(model_raw, "model", "crop_to_observed_range", True),
        interpolate=_require_bool(model_raw, "model", "interpolate", True),
        normalize_npol=_optional_integer_like(model_raw, "model", "normalize_npol"),
        normalize_parts=_optional_string(model_raw, "model", "normalize_parts"),
        instrument_resolution=_optional_number(model_raw, "model", "instrument_resolution"),
        rotational_broadening=_require_bool(model_raw, "model", "rotational_broadening", False),
        rotational_velocity=_optional_number(model_raw, "model", "rotational_velocity") or 0.0,
        doppler_shift=_require_bool(model_raw, "model", "doppler_shift", False),
        radial_velocity=_optional_number(model_raw, "model", "radial_velocity") or 0.0,
    )

    priors_raw = _require_mapping(root["priors"], "priors")
    priors: dict[str, PriorConfig] = {}
    for name, prior_raw in priors_raw.items():
        prior_section = f"priors.{name}"
        prior_data = _require_mapping(prior_raw, prior_section)
        _require_known_keys(prior_data, prior_section, PRIOR_KEYS)
        _require_keys(prior_data, prior_section, ("type",))
        priors[name] = PriorConfig(
            type=_require_string(prior_data, prior_section, "type"),
            value=_optional_number(prior_data, prior_section, "value"),
            min=_optional_number(prior_data, prior_section, "min"),
            max=_optional_number(prior_data, prior_section, "max"),
            mean=_optional_number(prior_data, prior_section, "mean"),
            sigma=_optional_number(prior_data, prior_section, "sigma"),
        )
    if not priors:
        raise ValueError("At least one prior must be configured")

    outputs_raw = _require_mapping(root["outputs"], "outputs")
    _require_known_keys(outputs_raw, "outputs", OUTPUT_KEYS)
    _require_keys(outputs_raw, "outputs", ("directory",))
    filenames_raw = outputs_raw.get("filenames", {})
    filenames = _require_mapping(filenames_raw, "outputs.filenames")
    for key, value in filenames.items():
        if not isinstance(value, str):
            raise ValueError(f"outputs.filenames.{key} must be a string.")
    outputs = OutputConfig(
        directory=_require_path(outputs_raw, "outputs", "directory", config_dir=config_dir),
        wavelength_type=_optional_string(outputs_raw, "outputs", "wavelength_type"),
        write_hdf5=_require_bool(outputs_raw, "outputs", "write_hdf5", True),
        write_plot=_require_bool(outputs_raw, "outputs", "write_plot", True),
        write_corner_plot=_require_bool(outputs_raw, "outputs", "write_corner_plot", True),
        write_bestfit_spectrum=_require_bool(outputs_raw, "outputs", "write_bestfit_spectrum", True),
        write_masked_observed=_require_bool(outputs_raw, "outputs", "write_masked_observed", True),
        write_run_summary=_require_bool(outputs_raw, "outputs", "write_run_summary", True),
        filenames=dict(filenames),
    )

    return AnalysisConfig(
        run=run,
        input_spectrum=input_spectrum,
        segments=segments,
        masking=masking,
        model=model,
        priors=priors,
        outputs=outputs,
        source_path=config_path,
    )


def load_config(path: str | Path) -> AnalysisConfig:
    config_path = absolute_user_path(path)
    raw = _load_json_with_duplicate_detection(config_path)
    return load_config_data(raw, config_path=config_path)
