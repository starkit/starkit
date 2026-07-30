from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from starkit.analysis.validation import verify_config


def _write_runtime_grid(path: Path) -> None:
    meta = pd.Series(
        {
            "parameters": json.dumps(["teff", "mh"], sort_keys=True),
            "wavelength_unit": json.dumps("Angstrom", sort_keys=True),
            "wavelength_type": json.dumps("vacuum", sort_keys=True),
            "flux_unit": json.dumps("erg / (cm2 s Angstrom)", sort_keys=True),
            "R": json.dumps(20000.0, sort_keys=True),
            "R_sampling": json.dumps(4.0, sort_keys=True),
            "format_version": json.dumps("v1.0", sort_keys=True),
            "grid_type": json.dumps("log", sort_keys=True),
        }
    )
    index = pd.DataFrame(
        {
            "teff": [5000.0, 5000.0, 5500.0, 5500.0],
            "mh": [0.0, 0.5, 0.0, 0.5],
        }
    )
    wavelength = pd.DataFrame({0: np.geomspace(1000.0, 1010.0, num=5)})

    meta.to_hdf(path, key="meta")
    index.to_hdf(path, key="index")
    wavelength.to_hdf(path, key="wavelength")
    with h5py.File(path, "a") as handle:
        handle.create_dataset(
            "fluxes",
            data=[
                [1.0, 1.1, 1.2, 1.3, 1.4],
                [0.95, 1.05, 1.15, 1.25, 1.35],
                [1.1, 1.2, 1.3, 1.4, 1.5],
                [1.05, 1.15, 1.25, 1.35, 1.45],
            ],
        )


def _write_spectrum(path: Path) -> None:
    data = np.column_stack(
        [
            np.geomspace(1000.0, 1010.0, num=5),
            np.array([1.0, 1.0, 0.98, 1.02, 1.01]),
            np.full(5, 0.05),
        ]
    )
    np.savetxt(path, data)


def _write_segments(path: Path, rows: list[tuple[float, float]]) -> None:
    np.savetxt(path, np.asarray(rows, dtype=float))


def _write_config(
    path: Path,
    *,
    grid_path: Path,
    spectrum_path: Path,
    priors: dict[str, dict[str, float | str]],
    model_overrides: dict[str, object] | None = None,
    segments_path: Path | None = None,
) -> None:
    model = {
        "likelihood": "chi2",
        "instrument_resolution": None,
        "rotational_broadening": False,
        "rotational_velocity": 0.0,
        "doppler_shift": False,
        "radial_velocity": 0.0,
        "crop_to_observed_range": True,
        "interpolate": True,
        "normalize_npol": None,
        "normalize_parts": None,
    }
    if model_overrides:
        model.update(model_overrides)

    config = {
        "run": {
            "name": "toy",
            "grid_path": str(grid_path),
            "grid_wavelength_type": "vacuum",
            "mpi": False,
            "mpi_processes": 1,
        },
        "input_spectrum": {
            "path": str(spectrum_path),
            "format": "ascii_columns",
            "wavelength_unit": "Angstrom",
            "flux_unit": "erg / (cm2 s Angstrom)",
            "wavelength_type": "vacuum",
            "columns": {
                "wavelength": 0,
                "flux": 1,
                "uncertainty": 2,
            },
            "skiprows": 0,
        },
        "segments": {
            "path": None if segments_path is None else str(segments_path),
            "mode": "file",
            "wavelength_reference_path": None,
            "wavelength_type": None,
            "resample_to_grid": False,
        },
        "model": model,
        "priors": priors,
        "outputs": {
            "directory": str(path.parent / "out"),
        },
    }
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")


def _codes(report) -> set[str]:
    return {issue.code for issue in report.issues}


def test_verify_config_reports_json_syntax_errors(tmp_path: Path) -> None:
    config_path = tmp_path / "bad.json"
    config_path.write_text('{"run": {"name": "toy",}', encoding="ascii")

    report = verify_config(config_path, deep=False)

    assert not report.ok
    assert "json.syntax" in _codes(report)


def test_verify_config_rejects_vrot_prior_without_rotational_broadening(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
            "vrot": {"type": "uniform", "min": 0.0, "max": 40.0},
        },
    )

    report = verify_config(config_path, deep=False)

    assert not report.ok
    assert "model.vrot.disabled" in _codes(report)


def test_verify_config_warns_when_vrad_is_left_fixed(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
        },
        model_overrides={
            "doppler_shift": True,
            "radial_velocity": 12.0,
        },
    )

    report = verify_config(config_path, deep=False)

    assert report.ok
    assert "model.vrad.fixed" in _codes(report)


def test_verify_config_rejects_r_prior_above_grid_resolution(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
            "R": {"type": "uniform", "min": 15000.0, "max": 20000.0},
        },
        model_overrides={
            "instrument_resolution": 18000.0,
        },
    )

    report = verify_config(config_path, deep=False)

    assert not report.ok
    assert "prior.R.above_grid_resolution" in _codes(report)


def test_verify_config_deep_check_catches_unknown_prior_name(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
            "bogus": {"type": "fixed", "value": 0.0},
        },
    )

    report = verify_config(config_path, deep=True)

    assert not report.ok
    assert report.deep_checks_run is True
    assert "deep.model" in _codes(report)
    assert any("Configured prior 'bogus'" in issue.message for issue in report.errors)


def test_verify_config_rejects_invalid_output_wavelength_type(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        priors={"teff": {"type": "uniform", "min": 5000.0, "max": 5500.0}},
    )
    config = json.loads(config_path.read_text(encoding="ascii"))
    config["outputs"]["wavelength_type"] = "water"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")

    report = verify_config(config_path, deep=False)

    assert not report.ok
    assert "outputs.wavelength_type" in _codes(report)


def test_verify_config_missing_path_mentions_config_directory(tmp_path: Path) -> None:
    config_dir = tmp_path / "cfg"
    config_dir.mkdir()
    inputs_dir = config_dir / "inputs"
    inputs_dir.mkdir()
    spectrum_path = inputs_dir / "spectrum.dat"
    _write_spectrum(spectrum_path)

    config_path = config_dir / "config.json"
    config = {
        "run": {
            "name": "toy",
            "grid_path": "inputs/missing_grid.h5",
            "mpi": False,
            "mpi_processes": 1,
        },
        "input_spectrum": {
            "path": "inputs/spectrum.dat",
            "format": "ascii_columns",
            "wavelength_unit": "Angstrom",
            "flux_unit": "erg / (cm2 s Angstrom)",
            "columns": {
                "wavelength": 0,
                "flux": 1,
            },
        },
        "model": {},
        "priors": {
            "teff": {
                "type": "uniform",
                "min": 5000.0,
                "max": 5500.0,
            }
        },
        "outputs": {
            "directory": "out",
        },
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")

    report = verify_config(config_path, deep=False)

    issue = next(issue for issue in report.errors if issue.code == "run.grid_path.missing")
    assert str(config_path.resolve()) in issue.message
    assert str(config_dir.resolve()) in issue.message
    assert "not the current working directory" in issue.message


def test_verify_config_rejects_normalize_parts_without_normalize_npol(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    segments_path = tmp_path / "segments.txt"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_segments(segments_path, [(1000.0, 1010.0)])
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        segments_path=segments_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
        },
        model_overrides={
            "normalize_parts": "segments",
        },
    )

    report = verify_config(config_path, deep=False)

    assert not report.ok
    assert "model.normalize_parts.requires_npol" in _codes(report)


def test_verify_config_rejects_normalize_parts_without_segments(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
        },
        model_overrides={
            "normalize_npol": 2,
            "normalize_parts": "segments",
        },
    )

    report = verify_config(config_path, deep=False)

    assert not report.ok
    assert "model.normalize_parts.requires_segments" in _codes(report)


def test_verify_config_allows_autodetect_without_segments(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
        },
        model_overrides={
            "normalize_npol": 2,
            "normalize_parts": "autodetect",
        },
    )

    report = verify_config(config_path, deep=False)

    assert report.ok
    assert "model.normalize_parts.requires_segments" not in _codes(report)


def test_verify_config_deep_check_rejects_too_short_normalize_parts_segment(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    segments_path = tmp_path / "segments.txt"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    _write_spectrum(spectrum_path)
    _write_segments(segments_path, [(1000.0, 1000.0), (1002.0, 1002.0)])
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        segments_path=segments_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
        },
        model_overrides={
            "normalize_npol": 2,
            "normalize_parts": "segments",
        },
    )

    report = verify_config(config_path, deep=True)

    assert not report.ok
    assert report.deep_checks_run is True
    assert "deep.normalize_parts" in _codes(report)


def test_verify_config_deep_check_rejects_too_short_autodetected_part(tmp_path: Path) -> None:
    grid_path = tmp_path / "grid.h5"
    spectrum_path = tmp_path / "spectrum.dat"
    config_path = tmp_path / "config.json"
    _write_runtime_grid(grid_path)
    np.savetxt(
        spectrum_path,
        np.column_stack(
            [
                np.array([1000.0, 1001.0, 1002.0, 2000.0]),
                np.array([1.0, 1.0, 0.98, 1.02]),
                np.full(4, 0.05),
            ]
        ),
    )
    _write_config(
        config_path,
        grid_path=grid_path,
        spectrum_path=spectrum_path,
        priors={
            "teff": {"type": "uniform", "min": 5000.0, "max": 5500.0},
        },
        model_overrides={
            "normalize_npol": 1,
            "normalize_parts": "autodetect",
        },
    )

    report = verify_config(config_path, deep=True)

    assert not report.ok
    assert report.deep_checks_run is True
    assert "deep.normalize_parts" in _codes(report)


def test_verify_config_uses_symlink_directory_for_relative_paths(tmp_path: Path) -> None:
    shared_dir = tmp_path / "shared"
    shared_dir.mkdir()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    config_path = shared_dir / "config.json"
    symlink_path = run_dir / "config.json"

    config = {
        "run": {
            "name": "toy",
            "grid_path": "inputs/missing_grid.h5",
            "mpi": False,
            "mpi_processes": 1,
        },
        "input_spectrum": {
            "path": "inputs/missing_spectrum.dat",
            "format": "ascii_columns",
            "wavelength_unit": "Angstrom",
            "flux_unit": "erg / (cm2 s Angstrom)",
            "columns": {
                "wavelength": 0,
                "flux": 1,
            },
        },
        "model": {},
        "priors": {
            "teff": {
                "type": "uniform",
                "min": 5000.0,
                "max": 5500.0,
            }
        },
        "outputs": {
            "directory": "out",
        },
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")
    symlink_path.symlink_to(config_path)

    report = verify_config(symlink_path, deep=False)

    grid_issue = next(issue for issue in report.errors if issue.code == "run.grid_path.missing")
    assert report.path == symlink_path
    assert str(run_dir / "inputs" / "missing_grid.h5") in grid_issue.message
    assert str(run_dir) in grid_issue.message
    assert str(shared_dir) not in grid_issue.message
