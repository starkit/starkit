from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import zipfile

import numpy as np
from astropy import units as u
from starkit.fix_spectrum1d import SKSpectrum1D
from starkit.utils.vacuumair_conversion import convert_air2vacuum, convert_vacuum2air

from .config import InputSpectrumConfig, MaskingConfig, SegmentsConfig


AUTODETECT_NORMALIZATION_GAP_FACTOR = 20.0


@dataclass(frozen=True)
class SpectrumArrays:
    wavelength: np.ndarray
    flux: np.ndarray
    uncertainty: np.ndarray | None
    include_regions: np.ndarray | None = None
    exclude_regions: np.ndarray | None = None
    segment_bounds: np.ndarray | None = None


def load_input_arrays(config: InputSpectrumConfig) -> SpectrumArrays:
    if config.format == "ascii_columns":
        return _load_ascii_columns(config)
    if config.format == "npy_columns":
        return _load_npy_columns(config)
    raise ValueError(f"Unsupported input format: {config.format}")


def _load_ascii_columns(config: InputSpectrumConfig) -> SpectrumArrays:
    data = np.loadtxt(config.path, skiprows=config.skiprows)
    if data.ndim != 2:
        raise ValueError(f"Expected 2D data in {config.path}")
    wavelength = data[:, _required_column(config, "wavelength")]
    flux = data[:, _required_column(config, "flux")]
    uncertainty = _load_uncertainty_from_columns_or_external(data, config)
    return SpectrumArrays(wavelength=wavelength, flux=flux, uncertainty=uncertainty)


def _load_npy_columns(config: InputSpectrumConfig) -> SpectrumArrays:
    data = np.load(config.path)
    if data.ndim != 2:
        raise ValueError(f"Expected 2D array in {config.path}")
    wavelength = data[:, _required_column(config, "wavelength")]
    flux = data[:, _required_column(config, "flux")]
    uncertainty = _load_uncertainty_from_columns_or_external(data, config)
    return SpectrumArrays(wavelength=wavelength, flux=flux, uncertainty=uncertainty)


def _required_column(config: InputSpectrumConfig, name: str) -> int:
    if name not in config.columns or config.columns[name] is None:
        raise ValueError(f"Column '{name}' must be configured for {config.format}")
    return int(config.columns[name])


def _load_uncertainty_from_columns_or_external(
    data: np.ndarray,
    config: InputSpectrumConfig,
) -> np.ndarray | None:
    uncertainty_idx = config.columns.get("uncertainty")
    if uncertainty_idx is not None:
        return data[:, int(uncertainty_idx)]

    if config.uncertainty_path is not None:
        uncertainty_data = _load_numeric_array(config.uncertainty_path)
        if config.uncertainty_column is None:
            if uncertainty_data.ndim == 1:
                return uncertainty_data
            raise ValueError("uncertainty_column is required for 2D uncertainty files")
        return uncertainty_data[:, int(config.uncertainty_column)]

    return None


def load_segments(
    config: SegmentsConfig,
    wavelength_unit: str | None = None,
    default_wavelength_type: str | None = None,
    target_wavelength_type: str | None = None,
) -> np.ndarray | None:
    if config.path is None:
        return None
    segments = np.loadtxt(config.path)
    if segments.ndim == 1:
        segments = segments.reshape(1, -1)
    if segments.shape[1] < 2:
        raise ValueError(f"Expected start/stop segment columns in {config.path}")
    return _convert_bounds(
        segments[:, :2],
        wavelength_unit=wavelength_unit,
        source_wavelength_type=_resolve_wavelength_type(config.wavelength_type, default_wavelength_type),
        target_wavelength_type=target_wavelength_type,
    )


def load_reference_wavelength(
    config: SegmentsConfig,
    wavelength_unit: str | None = None,
    default_wavelength_type: str | None = None,
    target_wavelength_type: str | None = None,
) -> np.ndarray | None:
    if config.wavelength_reference_path is None:
        return None
    data = _load_numeric_array(config.wavelength_reference_path)
    source_wavelength_type = _resolve_wavelength_type(config.wavelength_type, default_wavelength_type)
    if data.ndim == 1:
        return convert_wavelength_values(data, wavelength_unit, source_wavelength_type, target_wavelength_type)
    return convert_wavelength_values(
        data[:, 0],
        wavelength_unit,
        source_wavelength_type,
        target_wavelength_type,
    )


def _load_numeric_array(path: Path) -> np.ndarray:
    if _is_numpy_binary_file(path):
        loaded = np.load(path)
        if isinstance(loaded, np.lib.npyio.NpzFile):
            with loaded as archive:
                if len(archive.files) != 1:
                    raise ValueError(
                        f"Expected exactly one array in {path}, found {len(archive.files)}"
                    )
                return np.asarray(archive[archive.files[0]])
        return np.asarray(loaded)
    return np.asarray(np.loadtxt(path))


def _is_numpy_binary_file(path: Path) -> bool:
    if path.suffix.lower() in {".npy", ".npz"}:
        return True

    with path.open("rb") as handle:
        magic = handle.read(len(np.lib.format.MAGIC_PREFIX))
    if magic == np.lib.format.MAGIC_PREFIX:
        return True
    return zipfile.is_zipfile(path)


def build_analysis_grid(
    wavelength: np.ndarray,
    grid_resolution: float | None,
    grid_sampling: float | None,
    segments_config: SegmentsConfig,
    grid_reference_wavelength: np.ndarray | None = None,
    wavelength_unit: str | None = None,
    default_wavelength_type: str | None = None,
    target_wavelength_type: str | None = None,
    segments: np.ndarray | None = None,
) -> np.ndarray:
    if segments is None:
        segments = load_segments(
            segments_config,
            wavelength_unit=wavelength_unit,
            default_wavelength_type=default_wavelength_type,
            target_wavelength_type=target_wavelength_type,
        )
    if segments is None:
        return wavelength.copy()

    if not segments_config.resample_to_grid:
        masks = [(wavelength >= start) & (wavelength <= stop) for start, stop in segments]
        combined = np.logical_or.reduce(masks)
        return wavelength[combined]

    if grid_resolution is None or grid_sampling is None:
        raise ValueError("Grid R and R_sampling are required for resample_to_grid=true")

    reference = load_reference_wavelength(
        segments_config,
        wavelength_unit=wavelength_unit,
        default_wavelength_type=default_wavelength_type,
        target_wavelength_type=target_wavelength_type,
    )
    if reference is None:
        reference = grid_reference_wavelength
    if reference is None:
        reference = wavelength

    output_segments: list[np.ndarray] = []
    for start, stop in segments:
        indices = np.where((reference >= start) & (reference <= stop))[0]
        if indices.size == 0:
            continue
        seg_start = reference[indices[0]]
        seg_stop = reference[indices[-1]]
        segment_wave = np.exp(
            np.arange(
                np.log(seg_start),
                np.log(seg_stop),
                1.0 / (float(grid_sampling) * float(grid_resolution)),
            )
        )
        output_segments.append(segment_wave)

    if not output_segments:
        raise ValueError("No valid analysis segments intersected the wavelength grid")
    return np.concatenate(output_segments)


def prepare_analysis_arrays(
    raw_arrays: SpectrumArrays,
    input_config: InputSpectrumConfig,
    segments_config: SegmentsConfig,
    grid_resolution: float | None,
    grid_sampling: float | None,
    grid_reference_wavelength: np.ndarray | None = None,
    grid_wavelength_type: str | None = None,
) -> SpectrumArrays:
    analysis_wavelength_type = grid_wavelength_type or input_config.wavelength_type
    analysis_arrays = convert_wavelengths_to_grid_type(raw_arrays, input_config, grid_wavelength_type)
    segments = load_segments(
        segments_config,
        wavelength_unit=input_config.wavelength_unit,
        default_wavelength_type=input_config.wavelength_type,
        target_wavelength_type=analysis_wavelength_type,
    )
    analysis_wavelength = build_analysis_grid(
        wavelength=analysis_arrays.wavelength,
        grid_resolution=grid_resolution,
        grid_sampling=grid_sampling,
        segments_config=segments_config,
        grid_reference_wavelength=grid_reference_wavelength,
        wavelength_unit=input_config.wavelength_unit,
        default_wavelength_type=input_config.wavelength_type,
        target_wavelength_type=analysis_wavelength_type,
        segments=segments,
    )
    analysis_flux = resample_flux(
        analysis_arrays.wavelength,
        analysis_arrays.flux,
        analysis_wavelength,
    )
    analysis_uncertainty = None
    if analysis_arrays.uncertainty is not None:
        analysis_uncertainty = resample_flux(
            analysis_arrays.wavelength,
            analysis_arrays.uncertainty,
            analysis_wavelength,
        )

    return SpectrumArrays(
        wavelength=analysis_wavelength,
        flux=analysis_flux,
        uncertainty=analysis_uncertainty,
        segment_bounds=_merge_regions([segments]) if segments is not None else None,
    )


def resample_flux(
    source_wavelength: np.ndarray,
    source_flux: np.ndarray,
    target_wavelength: np.ndarray,
) -> np.ndarray:
    if target_wavelength.size == 0:
        return np.array([], dtype=float)

    source_min = float(np.min(source_wavelength))
    source_max = float(np.max(source_wavelength))
    target_min = float(np.min(target_wavelength))
    target_max = float(np.max(target_wavelength))
    tolerance = max(
        1e-8,
        1e-10 * max(abs(source_min), abs(source_max), abs(target_min), abs(target_max), 1.0),
    )
    if target_min < (source_min - tolerance) or target_max > (source_max + tolerance):
        raise ValueError(
            "Target wavelength grid extends outside the source wavelength coverage "
            f"({target_min} - {target_max} vs. {source_min} - {source_max})."
        )

    return np.interp(target_wavelength, source_wavelength, source_flux)


def apply_masking(
    wavelength: np.ndarray,
    flux: np.ndarray,
    uncertainty: np.ndarray | None,
    config: MaskingConfig,
    wavelength_unit: str | None = None,
    default_wavelength_type: str | None = None,
    target_wavelength_type: str | None = None,
    segment_bounds: np.ndarray | None = None,
) -> SpectrumArrays:
    source_wavelength_type = _resolve_wavelength_type(config.wavelength_type, default_wavelength_type)
    include_region_sets: list[np.ndarray] = []
    exclude_region_sets: list[np.ndarray] = []
    if config.include_regions_path is not None:
        include_region_sets.append(
            _load_region_file(
                config.include_regions_path,
                wavelength_unit=wavelength_unit,
                source_wavelength_type=source_wavelength_type,
                target_wavelength_type=target_wavelength_type,
            )
        )
    if config.exclude_regions_path is not None:
        exclude_region_sets.append(
            _load_region_file(
                config.exclude_regions_path,
                wavelength_unit=wavelength_unit,
                source_wavelength_type=source_wavelength_type,
                target_wavelength_type=target_wavelength_type,
            )
        )
    if config.regions_json_path is not None:
        json_regions = _load_regions_json(
            config.regions_json_path,
            wavelength_unit=wavelength_unit,
            source_wavelength_type=source_wavelength_type,
            target_wavelength_type=target_wavelength_type,
        )
        if len(json_regions["inclusion"]) > 0:
            include_region_sets.append(json_regions["inclusion"])
        if len(json_regions["exclusion"]) > 0:
            exclude_region_sets.append(json_regions["exclusion"])

    include_regions = _merge_regions(include_region_sets)
    exclude_regions = _merge_regions(exclude_region_sets)
    mask = np.ones_like(wavelength, dtype=bool)
    if include_regions.size > 0:
        mask = np.zeros_like(wavelength, dtype=bool)
        for start, stop in include_regions:
            mask |= (wavelength >= start) & (wavelength <= stop)

    for start, stop in exclude_regions:
        mask &= ~((wavelength >= start) & (wavelength <= stop))

    masked_uncertainty = None if uncertainty is None else uncertainty[mask]
    return SpectrumArrays(
        wavelength=wavelength[mask],
        flux=flux[mask],
        uncertainty=masked_uncertainty,
        include_regions=include_regions,
        exclude_regions=exclude_regions,
        segment_bounds=segment_bounds,
    )


def build_observed_spectrum(arrays: SpectrumArrays, config: InputSpectrumConfig) -> SKSpectrum1D:
    wavelength = arrays.wavelength * u.Unit(config.wavelength_unit)
    flux = arrays.flux * u.Unit(config.flux_unit)

    if arrays.uncertainty is not None:
        uncertainty = arrays.uncertainty * u.Unit(config.flux_unit)
    elif config.snr is not None:
        uncertainty = (np.full_like(arrays.flux, 1.0 / float(config.snr))) * flux.unit
    else:
        uncertainty = None

    spectrum = SKSpectrum1D.from_array(wavelength=wavelength, flux=flux, uncertainty=uncertainty)
    return spectrum.get_nan_cleaned()


def autodetect_normalization_segment_bounds(
    wavelength: np.ndarray,
    *,
    gap_factor: float = AUTODETECT_NORMALIZATION_GAP_FACTOR,
) -> np.ndarray:
    values = np.asarray(wavelength, dtype=float).reshape(-1)
    if values.size == 0:
        raise ValueError("Cannot autodetect normalization parts from an empty wavelength array.")
    if values.size == 1:
        return np.array([[values[0], values[0]]], dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("Cannot autodetect normalization parts from non-finite wavelengths.")
    if not np.all(values > 0.0):
        raise ValueError("Cannot autodetect normalization parts from non-positive wavelengths.")

    log_steps = np.diff(np.log(values))
    if not np.all(np.isfinite(log_steps)) or np.any(log_steps <= 0.0):
        raise ValueError("Cannot autodetect normalization parts from a wavelength grid that is not strictly increasing.")

    median_step = float(np.median(log_steps))
    if not np.isfinite(median_step) or median_step <= 0.0:
        raise ValueError("Cannot autodetect normalization parts because the wavelength step is invalid.")

    gap_indices = np.where(log_steps > float(gap_factor) * median_step)[0]
    part_starts = np.concatenate(([0], gap_indices + 1))
    part_stops = np.concatenate((gap_indices, [values.size - 1]))
    bounds = np.asarray(
        [[values[start], values[stop]] for start, stop in zip(part_starts, part_stops)],
        dtype=float,
    )
    return bounds.reshape(-1, 2)


def build_segment_normalization_parts(
    observed: SKSpectrum1D,
    segment_bounds: np.ndarray | None,
    npol: int,
) -> list[np.ndarray]:
    if segment_bounds is None:
        raise ValueError("Segment-based normalization requires analysis segments, but no segment bounds were provided.")

    bounds = np.asarray(segment_bounds, dtype=float).reshape(-1, 2)
    if bounds.size == 0:
        raise ValueError("Segment-based normalization requires at least one effective analysis segment.")

    wavelength = np.asarray(observed.wavelength.value, dtype=float)
    minimum_sample_count = int(npol) + 1
    parts: list[np.ndarray] = []
    for index, (start, stop) in enumerate(bounds, start=1):
        mask = (wavelength >= float(start)) & (wavelength <= float(stop))
        sample_count = int(mask.sum())
        if sample_count == 0:
            raise ValueError(
                f"Normalization segment {index} ({start:g} - {stop:g}) contains no observed samples after masking."
            )
        if sample_count < minimum_sample_count:
            raise ValueError(
                f"Normalization segment {index} ({start:g} - {stop:g}) contains {sample_count} sample(s), "
                f"but polynomial degree {npol} requires at least {minimum_sample_count}."
            )
        parts.append(mask)
    return parts


def _load_region_file(
    path: Path,
    wavelength_unit: str | None = None,
    source_wavelength_type: str | None = None,
    target_wavelength_type: str | None = None,
) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Region file does not exist: {path}")
    regions = np.loadtxt(path)
    if regions.ndim == 1:
        regions = regions.reshape(1, -1)
    if regions.shape[1] < 2:
        raise ValueError(f"Expected start/stop columns in region file {path}")
    return _convert_bounds(
        regions[:, :2],
        wavelength_unit=wavelength_unit,
        source_wavelength_type=source_wavelength_type,
        target_wavelength_type=target_wavelength_type,
    )


def _load_regions_json(
    path: Path,
    wavelength_unit: str | None = None,
    source_wavelength_type: str | None = None,
    target_wavelength_type: str | None = None,
) -> dict[str, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"Region JSON file does not exist: {path}")

    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)

    if isinstance(raw, list):
        raw_regions = raw
    elif isinstance(raw, dict):
        raw_regions = raw.get("regions")
    else:
        raise ValueError(f"Expected a JSON object or array in region JSON file {path}")

    if not isinstance(raw_regions, list):
        raise ValueError(f"Region JSON file {path} must contain a top-level 'regions' array")

    grouped = {"inclusion": [], "exclusion": []}
    for index, region in enumerate(raw_regions):
        if not isinstance(region, dict):
            raise ValueError(f"Region JSON entry {index + 1} in {path} is not an object")

        region_type = region.get("type")
        if region_type not in grouped:
            raise ValueError(
                f"Region JSON entry {index + 1} in {path} has invalid type {region_type!r}; "
                "expected 'inclusion' or 'exclusion'"
            )

        if "start" not in region or "end" not in region:
            raise ValueError(f"Region JSON entry {index + 1} in {path} is missing 'start' or 'end'")

        start = float(region["start"])
        end = float(region["end"])
        if not np.isfinite(start) or not np.isfinite(end):
            raise ValueError(f"Region JSON entry {index + 1} in {path} has non-finite bounds")

        grouped[region_type].append(_normalize_region_bounds(start, end))

    return {
        region_type: _convert_bounds(
            np.asarray(entries, dtype=float).reshape(-1, 2),
            wavelength_unit=wavelength_unit,
            source_wavelength_type=source_wavelength_type,
            target_wavelength_type=target_wavelength_type,
        )
        for region_type, entries in grouped.items()
    }


def _normalize_region_bounds(start: float, stop: float) -> tuple[float, float]:
    lower = float(min(start, stop))
    upper = float(max(start, stop))
    return lower, upper


def _merge_regions(region_sets: list[np.ndarray]) -> np.ndarray:
    if not region_sets:
        return np.empty((0, 2), dtype=float)

    merged_input = np.vstack([np.asarray(region_set, dtype=float).reshape(-1, 2) for region_set in region_sets])
    if merged_input.size == 0:
        return np.empty((0, 2), dtype=float)

    order = np.argsort(merged_input[:, 0], kind="mergesort")
    sorted_regions = merged_input[order]
    merged: list[list[float]] = []
    for start, stop in sorted_regions:
        if not merged or start > merged[-1][1]:
            merged.append([float(start), float(stop)])
            continue
        merged[-1][1] = max(merged[-1][1], float(stop))
    return np.asarray(merged, dtype=float).reshape(-1, 2)


def _convert_bounds(
    bounds: np.ndarray,
    wavelength_unit: str | None,
    source_wavelength_type: str | None,
    target_wavelength_type: str | None,
) -> np.ndarray:
    normalized_bounds = np.asarray(bounds, dtype=float).reshape(-1, 2)
    if normalized_bounds.size == 0:
        return normalized_bounds
    converted_bounds = convert_wavelength_values(
        normalized_bounds,
        wavelength_unit,
        source_wavelength_type,
        target_wavelength_type,
    )
    return np.asarray(
        [_normalize_region_bounds(start, stop) for start, stop in converted_bounds],
        dtype=float,
    ).reshape(-1, 2)


def _resolve_wavelength_type(configured_type: str | None, default_type: str | None) -> str | None:
    if configured_type is not None:
        return configured_type
    return default_type


def convert_wavelengths_to_grid_type(
    arrays: SpectrumArrays,
    config: InputSpectrumConfig,
    grid_wavelength_type: str | None,
) -> SpectrumArrays:
    return SpectrumArrays(
        wavelength=convert_wavelength_values(
            arrays.wavelength,
            config.wavelength_unit,
            config.wavelength_type,
            grid_wavelength_type,
        ),
        flux=arrays.flux,
        uncertainty=arrays.uncertainty,
        include_regions=arrays.include_regions,
        exclude_regions=arrays.exclude_regions,
        segment_bounds=arrays.segment_bounds,
    )


def convert_wavelength_values(
    values: np.ndarray,
    wavelength_unit: str | None,
    source_type: str | None,
    target_type: str | None,
) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if source_type is None or target_type is None or source_type == target_type:
        return array

    if wavelength_unit is None:
        raise ValueError("wavelength_unit is required when converting wavelength values")
    if source_type not in {"air", "vacuum"}:
        raise ValueError(f"Unsupported source wavelength_type: {source_type}")
    if target_type not in {"air", "vacuum"}:
        raise ValueError(f"Unsupported target wavelength_type: {target_type}")

    wavelength_quantity = array * u.Unit(wavelength_unit)
    if source_type == "air" and target_type == "vacuum":
        converted = convert_air2vacuum(wavelength_quantity)
    elif source_type == "vacuum" and target_type == "air":
        converted = convert_vacuum2air(wavelength_quantity)
    else:
        raise ValueError(f"Unsupported wavelength conversion: {source_type} -> {target_type}")
    return converted.to_value(u.Unit(wavelength_unit))
