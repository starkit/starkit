from __future__ import annotations

import json
from pathlib import Path
import warnings

import h5py
import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning

from starkit.gridkit.inspector import format_inspection_summary, inspect_grid


def _write_test_grid(
    path: Path,
    *,
    legacy_meta: bool = False,
    wavelength_unit: str = "Angstrom",
    flux_unit: str | None = None,
    grid_type: str | None = "log",
    log_spaced: bool = True,
) -> None:
    raw_meta = {
        "parameters": ["teff", "logg", "mh"],
        "wavelength_unit": wavelength_unit,
        "wavelength_type": "vacuum",
        "flux_unit": flux_unit or f"erg / (cm2 s {wavelength_unit})",
        "R": 20000.0,
        "R_sampling": 4.0,
        "format_version": [1, 0],
    }
    if grid_type is not None:
        raw_meta["grid_type"] = grid_type
    meta = pd.Series(
        raw_meta if legacy_meta else {key: json.dumps(value, sort_keys=True) for key, value in raw_meta.items()}
    )
    index = pd.DataFrame(
        {
            "teff": [5000.0, 5000.0, 5500.0],
            "logg": [4.0, 4.5, 4.0],
            "mh": [0.0, 0.0, 0.0],
        }
    )
    wavelength_array = np.geomspace(1000.0, 1010.0, num=3) if log_spaced else [1000.0, 1005.0, 1010.0]
    wavelength = pd.DataFrame({0: wavelength_array})
    if legacy_meta:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PerformanceWarning)
            meta.to_hdf(path, key="meta")
    else:
        meta.to_hdf(path, key="meta")
    index.to_hdf(path, key="index")
    wavelength.to_hdf(path, key="wavelength")
    with h5py.File(path, "a") as handle:
        handle.create_dataset(
            "fluxes",
            data=[
                [1.0, 1.1, 1.2],
                [0.9, 1.0, 1.1],
                [1.05, 1.15, 1.25],
            ],
        )


def test_inspect_grid_and_format_summary(tmp_path: Path) -> None:
    grid_path = tmp_path / "toy_grid.h5"
    _write_test_grid(grid_path)

    inspection = inspect_grid(grid_path)

    assert inspection.parameter_names == ["teff", "logg", "mh"]
    assert inspection.point_count == 3
    assert inspection.expected_cartesian_point_count == 4
    assert inspection.occupancy_fraction == 0.75
    assert inspection.wavelength_type == "vacuum"
    assert inspection.grid_type == "log"
    assert inspection.wavelength_spacing == "log"
    assert inspection.starkit_compatible is True
    assert inspection.resolution == 20000.0
    assert inspection.resolution_sampling == 4.0
    assert inspection.grid_shape == (3, 3)
    assert any("not a full cartesian product" in warning for warning in inspection.warnings)
    summary = format_inspection_summary(inspection)
    assert "- teff: 5000 to 5500 (2 unique; sample: 5000, 5500)" in summary
    assert "5.0000e+03" not in summary
    assert "Wavelength range: 1000.000 to 1010.000 Angstrom" in summary
    assert "1.0000e+03 to 1.0100e+03" not in summary
    assert "Grid type: log" in summary
    assert "Wavelength spacing: log" in summary
    assert "Starkit compatible: yes" in summary
    assert "R: 20000" in summary
    assert "2.0000e+04" not in summary

def test_inspect_grid_supports_legacy_meta_series(tmp_path: Path) -> None:
    grid_path = tmp_path / "legacy_grid.h5"
    _write_test_grid(grid_path, legacy_meta=True)

    inspection = inspect_grid(grid_path)

    assert inspection.parameter_names == ["teff", "logg", "mh"]
    assert inspection.wavelength_type == "vacuum"
    assert inspection.grid_type == "log"
    assert inspection.resolution == 20000.0


def test_inspection_summary_formats_nanometer_wavelength_range_without_scientific_notation(tmp_path: Path) -> None:
    grid_path = tmp_path / "toy_grid_nm.h5"
    _write_test_grid(grid_path, wavelength_unit="nm")

    inspection = inspect_grid(grid_path)
    summary = format_inspection_summary(inspection)

    assert "Wavelength range: 1000.000 to 1010.000 nm" in summary
    assert "1.0000e+03 to 1.0100e+03" not in summary


def test_inspection_summary_explains_dimensionless_flux_unit(tmp_path: Path) -> None:
    grid_path = tmp_path / "toy_grid_dimensionless.h5"
    _write_test_grid(grid_path, flux_unit="1")

    inspection = inspect_grid(grid_path)
    summary = format_inspection_summary(inspection)

    assert "Flux unit: 1 (dimensionless)" in summary


def test_inspect_grid_warns_for_missing_grid_type_and_linear_wavelength_spacing(tmp_path: Path) -> None:
    grid_path = tmp_path / "toy_grid_linear.h5"
    _write_test_grid(grid_path, grid_type=None, log_spaced=False)

    inspection = inspect_grid(grid_path)

    assert inspection.grid_type is None
    assert inspection.wavelength_spacing == "linear"
    assert inspection.starkit_compatible is False
    assert any("does not declare grid_type" in warning for warning in inspection.warnings)
    assert any("not evenly spaced in log(lambda)" in warning for warning in inspection.warnings)

def test_inspect_grid_warns_when_grid_type_claims_log_but_sampling_is_linear(tmp_path: Path) -> None:
    grid_path = tmp_path / "toy_grid_mismatched.h5"
    _write_test_grid(grid_path, grid_type="log", log_spaced=False)

    inspection = inspect_grid(grid_path)

    assert inspection.grid_type == "log"
    assert inspection.wavelength_spacing == "linear"
    assert inspection.starkit_compatible is False
    assert any("declares grid_type 'log', but the wavelength sampling looks linear" in warning for warning in inspection.warnings)
