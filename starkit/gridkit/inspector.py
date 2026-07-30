from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd


@dataclass(frozen=True)
class ParameterInspection:
    name: str
    minimum: float
    maximum: float
    unique_count: int
    unique_values_preview: list[float]


@dataclass(frozen=True)
class GridInspection:
    path: str
    parameter_names: list[str]
    parameters: list[ParameterInspection]
    point_count: int
    expected_cartesian_point_count: int
    occupancy_fraction: float | None
    wavelength_count: int
    wavelength_min: float | None
    wavelength_max: float | None
    wavelength_unit: str | None
    wavelength_type: str | None
    grid_type: str | None
    wavelength_spacing: str | None
    starkit_compatible: bool
    flux_unit: str | None
    resolution: float | None
    resolution_sampling: float | None
    grid_shape: tuple[int, ...] | None
    format_version: Any
    warnings: list[str]

    @property
    def flux_shape(self) -> tuple[int, ...] | None:
        return self.grid_shape


def inspect_grid(path: str | Path) -> GridInspection:
    grid_path = Path(path).expanduser().resolve()
    meta = pd.read_hdf(grid_path, "meta")
    index = pd.read_hdf(grid_path, "index")
    wavelength_frame = pd.read_hdf(grid_path, "wavelength")

    with h5py.File(grid_path, "r") as handle:
        flux_dataset = handle.get("fluxes")
        grid_shape = tuple(flux_dataset.shape) if flux_dataset is not None else None

    meta_map = _series_to_dict(meta)
    parameter_names = _normalize_parameters(meta_map.get("parameters", []))

    warnings: list[str] = []
    parameters: list[ParameterInspection] = []
    unique_count_product = 1

    for parameter_name in parameter_names:
        if parameter_name not in index.columns:
            warnings.append(
                f"Grid metadata declares parameter '{parameter_name}', but the index table does not contain that column."
            )
            continue

        unique_values = np.sort(index[parameter_name].dropna().unique())
        if unique_values.size == 0:
            warnings.append(f"Parameter '{parameter_name}' has no non-null sampled values in the index.")
            continue

        unique_count_product *= int(unique_values.size)
        minimum = float(unique_values[0])
        maximum = float(unique_values[-1])
        if unique_values.size == 1:
            warnings.append(
                "Parameter "
                f"'{parameter_name}' has only one sampled value ({_format_parameter_value(parameter_name, minimum)}); "
                "its prior will be fixed."
            )

        parameters.append(
            ParameterInspection(
                name=parameter_name,
                minimum=minimum,
                maximum=maximum,
                unique_count=int(unique_values.size),
                unique_values_preview=[float(value) for value in unique_values[:10]],
            )
        )

    point_count = int(len(index))
    expected_cartesian_point_count = unique_count_product if parameters else 0
    occupancy_fraction = _compute_occupancy_fraction(point_count, expected_cartesian_point_count)

    if occupancy_fraction is not None and occupancy_fraction < 1.0:
        warnings.append(
            "Grid index is not a full cartesian product of the sampled parameter values "
            f"(occupancy {_format_fraction(occupancy_fraction)})."
        )
    if occupancy_fraction is not None and occupancy_fraction < 0.7:
        warnings.append(
            "Grid occupancy is substantially below a full cartesian product; interpolation coverage may be irregular."
        )

    wavelength_values = _extract_wavelength_values(wavelength_frame)
    wavelength_min = float(np.min(wavelength_values)) if wavelength_values.size else None
    wavelength_max = float(np.max(wavelength_values)) if wavelength_values.size else None
    wavelength_count = int(wavelength_values.size)

    wavelength_type = _normalize_optional_string(meta_map.get("wavelength_type"))
    if wavelength_type is None:
        warnings.append("Grid metadata does not declare a wavelength_type; downstream tools may need an explicit override.")
    elif wavelength_type not in {"air", "vacuum"}:
        warnings.append(f"Grid wavelength_type '{wavelength_type}' is not one of the expected values: air, vacuum.")

    resolution = _maybe_float(meta_map.get("R"))
    resolution_sampling = _maybe_float(meta_map.get("R_sampling"))
    if resolution is None:
        warnings.append("Grid metadata does not provide R; resolution-dependent defaults cannot be inferred.")
    if resolution_sampling is None:
        warnings.append("Grid metadata does not provide R_sampling; log-grid resampling defaults may need adjustment.")

    grid_type = _normalize_optional_string(meta_map.get("grid_type"))
    wavelength_spacing = _classify_wavelength_spacing(wavelength_values)
    if grid_type is None:
        warnings.append("Grid metadata does not declare grid_type; upstream Starkit requires grid_type='log'.")
    elif grid_type != "log":
        warnings.append(f"Grid metadata declares grid_type '{grid_type}', but upstream Starkit only supports 'log'.")

    if wavelength_spacing == "invalid":
        warnings.append("Grid wavelength array is not strictly increasing and positive; Starkit cannot treat it as a log grid.")
    elif wavelength_spacing != "log":
        warnings.append(
            "Grid wavelength array is not evenly spaced in log(lambda); upstream Starkit expects a log-spaced wavelength grid."
        )

    if grid_type == "log" and wavelength_spacing not in {None, "log"}:
        warnings.append(
            f"Grid metadata declares grid_type 'log', but the wavelength sampling looks {wavelength_spacing}."
        )

    starkit_compatible = bool(grid_type == "log" and wavelength_spacing == "log")

    if not parameter_names:
        warnings.append("Grid metadata does not list any interpolation parameters.")

    return GridInspection(
        path=str(grid_path),
        parameter_names=parameter_names,
        parameters=parameters,
        point_count=point_count,
        expected_cartesian_point_count=expected_cartesian_point_count,
        occupancy_fraction=occupancy_fraction,
        wavelength_count=wavelength_count,
        wavelength_min=wavelength_min,
        wavelength_max=wavelength_max,
        wavelength_unit=_normalize_optional_string(meta_map.get("wavelength_unit")),
        wavelength_type=wavelength_type,
        grid_type=grid_type,
        wavelength_spacing=wavelength_spacing,
        starkit_compatible=starkit_compatible,
        flux_unit=_normalize_optional_string(meta_map.get("flux_unit")),
        resolution=resolution,
        resolution_sampling=resolution_sampling,
        grid_shape=grid_shape,
        format_version=meta_map.get("format_version"),
        warnings=warnings,
    )


def format_inspection_summary(inspection: GridInspection) -> str:
    lines = [
        f"Grid: {inspection.path}",
        f"Points: {inspection.point_count}",
        f"Parameters: {', '.join(inspection.parameter_names) if inspection.parameter_names else 'none'}",
    ]

    if inspection.expected_cartesian_point_count:
        lines.append(
            "Occupancy: "
            f"{inspection.point_count}/{inspection.expected_cartesian_point_count}"
            + (
                f" ({_format_fraction(inspection.occupancy_fraction)})"
                if inspection.occupancy_fraction is not None
                else ""
            )
        )

    lines.extend(
        [
            f"Wavelengths: {inspection.wavelength_count}",
            "Wavelength range: "
            f"{_format_optional_wavelength(inspection.wavelength_min, inspection.wavelength_unit)} "
            f"to {_format_optional_wavelength(inspection.wavelength_max, inspection.wavelength_unit)} "
            f"{inspection.wavelength_unit or ''}".rstrip(),
            f"Wavelength type: {inspection.wavelength_type or 'unknown'}",
            f"Grid type: {inspection.grid_type or 'unknown'}",
            f"Wavelength spacing: {inspection.wavelength_spacing or 'unknown'}",
            f"Starkit compatible: {'yes' if inspection.starkit_compatible else 'no'}",
            f"Flux unit: {_format_flux_unit(inspection.flux_unit)}",
            f"R: {_format_optional_resolution(inspection.resolution)}",
            f"R_sampling: {_format_optional_float(inspection.resolution_sampling)}",
            f"Grid shape: {inspection.grid_shape or 'unknown'}",
        ]
    )

    if inspection.parameters:
        lines.append("Parameter sampling:")
        for parameter in inspection.parameters:
            preview = ", ".join(
                _format_parameter_value(parameter.name, value) for value in parameter.unique_values_preview
            )
            suffix = " ..." if parameter.unique_count > len(parameter.unique_values_preview) else ""
            lines.append(
                "  "
                f"- {parameter.name}: "
                f"{_format_parameter_value(parameter.name, parameter.minimum)} "
                f"to {_format_parameter_value(parameter.name, parameter.maximum)} "
                f"({parameter.unique_count} unique; sample: {preview}{suffix})"
            )

    if inspection.warnings:
        lines.append("Warnings:")
        for warning in inspection.warnings:
            lines.append(f"  - {warning}")

    return "\n".join(lines)


def inspection_to_dict(inspection: GridInspection) -> dict[str, Any]:
    return asdict(inspection)


def inspection_to_json(inspection: GridInspection) -> str:
    return json.dumps(inspection_to_dict(inspection), indent=2, sort_keys=True)


def _format_flux_unit(flux_unit: str | None) -> str:
    if flux_unit is None:
        return "unknown"
    if flux_unit == "1":
        return "1 (dimensionless)"
    return flux_unit


def _series_to_dict(meta: pd.Series | pd.DataFrame) -> dict[str, Any]:
    if isinstance(meta, pd.DataFrame):
        if meta.shape[1] == 1:
            return _decode_meta_mapping(meta.iloc[:, 0].to_dict())
        return _decode_meta_mapping(meta.iloc[0].to_dict())
    return _decode_meta_mapping(meta.to_dict())


def _decode_meta_mapping(raw_meta: dict[str, Any]) -> dict[str, Any]:
    return {key: _decode_meta_value(value) for key, value in raw_meta.items()}


def _decode_meta_value(value: Any) -> Any:
    if not isinstance(value, str):
        return value
    try:
        return json.loads(value)
    except json.JSONDecodeError:
        return value


def _normalize_parameters(raw_parameters: Any) -> list[str]:
    if raw_parameters is None:
        return []
    if isinstance(raw_parameters, str):
        return [raw_parameters]
    if isinstance(raw_parameters, (list, tuple)):
        return [str(item) for item in raw_parameters]
    if isinstance(raw_parameters, np.ndarray):
        return [str(item) for item in raw_parameters.tolist()]
    return [str(raw_parameters)]


def _extract_wavelength_values(frame: pd.DataFrame | pd.Series) -> np.ndarray:
    if isinstance(frame, pd.Series):
        values = frame.to_numpy()
    else:
        values = frame.to_numpy()
    return np.asarray(values).reshape(-1)


def _classify_wavelength_spacing(values: np.ndarray) -> str | None:
    wavelength_values = np.asarray(values, dtype=float).reshape(-1)
    if wavelength_values.size < 2:
        return None
    if not np.all(np.isfinite(wavelength_values)):
        return "invalid"

    diffs = np.diff(wavelength_values)
    if np.any(diffs <= 0.0) or np.any(wavelength_values <= 0.0):
        return "invalid"

    if _is_evenly_spaced(np.diff(np.log(wavelength_values))):
        return "log"
    if _is_evenly_spaced(diffs):
        return "linear"
    return "irregular"


def _is_evenly_spaced(differences: np.ndarray) -> bool:
    if differences.size == 0:
        return True
    reference = float(differences[0])
    tolerance = max(1e-12, abs(reference) * 1e-4)
    return bool(np.all(np.abs(differences - reference) <= tolerance))


def _compute_occupancy_fraction(point_count: int, expected_count: int) -> float | None:
    if expected_count <= 0:
        return None
    return float(point_count) / float(expected_count)


def _maybe_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and not value.strip():
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if np.isnan(result):
        return None
    return result


def _normalize_optional_string(value: Any) -> str | None:
    if value is None:
        return None
    result = str(value).strip()
    return result or None


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return "unknown"
    return _format_float(value)


def _format_optional_resolution(value: float | None) -> str:
    if value is None:
        return "unknown"
    return _format_resolution(value)


def _format_optional_wavelength(value: float | None, unit: str | None) -> str:
    if value is None:
        return "unknown"
    return _format_wavelength(value, unit)


def _format_parameter_value(parameter_name: str, value: float) -> str:
    if parameter_name == "teff":
        return _format_teff(value)
    return _format_float(value)


def _format_float(value: float) -> str:
    if not np.isfinite(value):
        return str(value)
    if value == 0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1_000 or magnitude < 0.01:
        return f"{value:.4e}"
    return f"{value:.6f}".rstrip("0").rstrip(".")


def _format_teff(value: float) -> str:
    if not np.isfinite(value):
        return str(value)
    rounded = round(value)
    if np.isclose(value, rounded):
        return str(int(rounded))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def _format_resolution(value: float) -> str:
    if not np.isfinite(value):
        return str(value)
    rounded = round(value)
    if np.isclose(value, rounded):
        return str(int(rounded))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def _format_wavelength(value: float, unit: str | None) -> str:
    if not np.isfinite(value):
        return str(value)
    if unit is not None and unit.strip().lower() in {"angstrom", "angstroms", "nanometer", "nanometers", "nm"}:
        return f"{value:.3f}"
    return _format_float(value)


def _format_fraction(value: float | None) -> str:
    if value is None:
        return "unknown"
    return f"{value:.3f}"
