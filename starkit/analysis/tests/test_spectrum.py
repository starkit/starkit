from __future__ import annotations

import json
from pathlib import Path

from astropy import units as u
import numpy as np
import pytest
from starkit.fix_spectrum1d import SKSpectrum1D
from starkit.utils.vacuumair_conversion import convert_air2vacuum

from starkit.analysis.config import InputSpectrumConfig, MaskingConfig, SegmentsConfig
from starkit.analysis.spectrum import (
    SpectrumArrays,
    apply_masking,
    autodetect_normalization_segment_bounds,
    build_segment_normalization_parts,
    load_input_arrays,
    load_reference_wavelength,
    prepare_analysis_arrays,
    resample_flux,
)


def test_resample_flux_rejects_targets_outside_source_range() -> None:
    with pytest.raises(ValueError, match="outside the source wavelength coverage"):
        resample_flux(
            np.array([100.0, 101.0]),
            np.array([1.0, 2.0]),
            np.array([99.5, 100.5, 101.5]),
        )


def test_apply_masking_combines_text_and_json_regions(tmp_path: Path) -> None:
    include_path = tmp_path / "include.txt"
    include_path.write_text("2 3\n", encoding="ascii")

    exclude_path = tmp_path / "exclude.txt"
    exclude_path.write_text("3 3\n", encoding="ascii")

    regions_json_path = tmp_path / "regions.json"
    regions_json_path.write_text(
        json.dumps(
            {
                "format": "selection_tool_regions",
                "version": 1,
                "sourceName": "ignored.dat",
                "regions": [
                    {"type": "inclusion", "start": 6, "end": 5},
                    {"type": "exclusion", "start": 6, "end": 6},
                ],
            }
        )
        + "\n",
        encoding="utf-8",
    )

    arrays = apply_masking(
        wavelength=np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]),
        flux=np.array([10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0]),
        uncertainty=np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]),
        config=MaskingConfig(
            include_regions_path=include_path,
            exclude_regions_path=exclude_path,
            regions_json_path=regions_json_path,
        ),
    )

    np.testing.assert_allclose(arrays.wavelength, np.array([2.0, 5.0]))
    np.testing.assert_allclose(arrays.flux, np.array([20.0, 50.0]))
    np.testing.assert_allclose(arrays.uncertainty, np.array([1.0, 1.0]))
    np.testing.assert_allclose(arrays.include_regions, np.array([[2.0, 3.0], [5.0, 6.0]]))
    np.testing.assert_allclose(arrays.exclude_regions, np.array([[3.0, 3.0], [6.0, 6.0]]))


def test_prepare_analysis_arrays_defaults_segment_wavelength_type_to_observed_spectrum(
    tmp_path: Path,
) -> None:
    segments_path = tmp_path / "segments.txt"
    segments_path.write_text("5000 5002\n", encoding="ascii")

    raw_arrays = SpectrumArrays(
        wavelength=np.array([5000.0, 5001.0, 5002.0]),
        flux=np.array([1.0, 0.9, 1.1]),
        uncertainty=np.array([0.1, 0.1, 0.1]),
    )
    input_config = InputSpectrumConfig(
        path=tmp_path / "spectrum.dat",
        format="ascii_columns",
        wavelength_unit="Angstrom",
        flux_unit="erg / (cm2 s Angstrom)",
        wavelength_type="air",
    )
    segments_config = SegmentsConfig(
        path=segments_path,
        resample_to_grid=False,
    )
    expected_wavelength = convert_air2vacuum(raw_arrays.wavelength * u.angstrom).to_value(u.angstrom)

    analysis = prepare_analysis_arrays(
        raw_arrays=raw_arrays,
        input_config=input_config,
        segments_config=segments_config,
        grid_resolution=5000.0,
        grid_sampling=1.0,
        grid_reference_wavelength=None,
        grid_wavelength_type="vacuum",
    )

    np.testing.assert_allclose(analysis.wavelength, expected_wavelength)
    np.testing.assert_allclose(analysis.flux, raw_arrays.flux)
    np.testing.assert_allclose(analysis.uncertainty, raw_arrays.uncertainty)
    np.testing.assert_allclose(
        analysis.segment_bounds,
        np.array([[expected_wavelength[0], expected_wavelength[-1]]]),
    )


def test_apply_masking_defaults_region_wavelength_type_to_observed_spectrum(tmp_path: Path) -> None:
    exclude_path = tmp_path / "exclude.txt"
    exclude_path.write_text("5001 5001\n", encoding="ascii")

    wavelength = convert_air2vacuum(np.array([5000.0, 5001.0, 5002.0]) * u.angstrom).to_value(u.angstrom)
    arrays = apply_masking(
        wavelength=wavelength,
        flux=np.array([10.0, 20.0, 30.0]),
        uncertainty=np.array([1.0, 1.0, 1.0]),
        config=MaskingConfig(exclude_regions_path=exclude_path),
        wavelength_unit="Angstrom",
        default_wavelength_type="air",
        target_wavelength_type="vacuum",
    )

    np.testing.assert_allclose(arrays.wavelength, np.array([wavelength[0], wavelength[2]]))
    np.testing.assert_allclose(arrays.flux, np.array([10.0, 30.0]))
    np.testing.assert_allclose(arrays.uncertainty, np.array([1.0, 1.0]))


def test_load_reference_wavelength_accepts_text_files(tmp_path: Path) -> None:
    reference_path = tmp_path / "reference.txt"
    reference_path.write_text("5000\n5001\n5002\n", encoding="ascii")

    loaded = load_reference_wavelength(
        SegmentsConfig(wavelength_reference_path=reference_path),
    )

    np.testing.assert_allclose(loaded, np.array([5000.0, 5001.0, 5002.0]))


def test_load_reference_wavelength_accepts_extensionless_npz_file(tmp_path: Path) -> None:
    reference_archive = tmp_path / "reference.npz"
    reference_path = tmp_path / "reference.dat"
    np.savez(reference_archive, wavelength=np.array([5000.0, 5001.0, 5002.0]))
    reference_archive.rename(reference_path)

    loaded = load_reference_wavelength(
        SegmentsConfig(wavelength_reference_path=reference_path),
    )

    np.testing.assert_allclose(loaded, np.array([5000.0, 5001.0, 5002.0]))


def test_load_input_arrays_accepts_extensionless_npy_uncertainty_file(tmp_path: Path) -> None:
    spectrum_path = tmp_path / "spectrum.dat"
    np.savetxt(
        spectrum_path,
        np.column_stack(
            [
                np.array([5000.0, 5001.0, 5002.0]),
                np.array([1.0, 0.9, 1.1]),
            ]
        ),
    )
    uncertainty_path = tmp_path / "uncertainty.dat"
    with uncertainty_path.open("wb") as handle:
        np.save(handle, np.array([0.1, 0.2, 0.3]))

    arrays = load_input_arrays(
        InputSpectrumConfig(
            path=spectrum_path,
            format="ascii_columns",
            wavelength_unit="Angstrom",
            flux_unit="erg / (cm2 s Angstrom)",
            columns={"wavelength": 0, "flux": 1},
            uncertainty_path=uncertainty_path,
        )
    )

    np.testing.assert_allclose(arrays.wavelength, np.array([5000.0, 5001.0, 5002.0]))
    np.testing.assert_allclose(arrays.flux, np.array([1.0, 0.9, 1.1]))
    np.testing.assert_allclose(arrays.uncertainty, np.array([0.1, 0.2, 0.3]))


def test_build_segment_normalization_parts_uses_effective_segment_bounds() -> None:
    observed = SKSpectrum1D.from_array(
        wavelength=np.array([5000.0, 5001.0, 5002.0, 5005.0, 5006.0]) * u.angstrom,
        flux=np.array([1.0, 0.9, 1.1, 1.0, 0.95]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=np.full(5, 0.1) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
    )

    parts = build_segment_normalization_parts(
        observed,
        np.array([[5000.0, 5002.0], [5005.0, 5006.0]]),
        npol=1,
    )

    assert len(parts) == 2
    np.testing.assert_array_equal(parts[0], np.array([True, True, True, False, False]))
    np.testing.assert_array_equal(parts[1], np.array([False, False, False, True, True]))


def test_build_segment_normalization_parts_rejects_short_segments() -> None:
    observed = SKSpectrum1D.from_array(
        wavelength=np.array([5000.0, 5001.0, 5005.0]) * u.angstrom,
        flux=np.array([1.0, 0.9, 1.1]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=np.full(3, 0.1) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
    )

    with pytest.raises(ValueError, match="requires at least 3"):
        build_segment_normalization_parts(
            observed,
            np.array([[5000.0, 5001.0], [5005.0, 5005.0]]),
            npol=2,
        )


def test_autodetect_normalization_segment_bounds_returns_single_part_without_large_gaps() -> None:
    bounds = autodetect_normalization_segment_bounds(
        np.array([5000.0, 5000.1, 5000.2, 5000.3]),
    )

    np.testing.assert_allclose(bounds, np.array([[5000.0, 5000.3]]))


def test_autodetect_normalization_segment_bounds_detects_large_gaps() -> None:
    bounds = autodetect_normalization_segment_bounds(
        np.array([5000.0, 5000.1, 5000.2, 5010.0, 5010.1, 5030.0]),
    )

    np.testing.assert_allclose(
        bounds,
        np.array([
            [5000.0, 5000.2],
            [5010.0, 5010.1],
            [5030.0, 5030.0],
        ]),
    )


def test_autodetect_bounds_do_not_split_on_masking_holes() -> None:
    bounds = autodetect_normalization_segment_bounds(
        np.array([5000.0, 5000.1, 5000.2, 5000.3, 5000.4]),
    )
    observed = SKSpectrum1D.from_array(
        wavelength=np.array([5000.0, 5000.1, 5000.4]) * u.angstrom,
        flux=np.array([1.0, 0.9, 1.1]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=np.full(3, 0.1) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
    )

    parts = build_segment_normalization_parts(observed, bounds, npol=1)

    assert len(parts) == 1
    np.testing.assert_array_equal(parts[0], np.array([True, True, True]))
