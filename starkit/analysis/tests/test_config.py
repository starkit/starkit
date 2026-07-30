from __future__ import annotations

import json
from pathlib import Path

import pytest

from starkit.analysis.config import load_config


def _base_config() -> dict[str, object]:
    return {
        "run": {
            "name": "toy",
            "grid_path": "grid.h5",
            "mpi": False,
            "mpi_processes": 1,
        },
        "input_spectrum": {
            "path": "spectrum.dat",
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
                "min": 3000.0,
                "max": 5000.0,
            }
        },
        "outputs": {
            "directory": "out",
        },
    }


def test_load_config_resolves_relative_paths_against_config_file(tmp_path: Path) -> None:
    config_dir = tmp_path / "config_dir"
    config_dir.mkdir()
    config_path = config_dir / "config.json"

    config = _base_config()
    config["run"] = {
        "name": "toy",
        "grid_path": "inputs/grid.h5",
        "mpi": False,
        "mpi_processes": 1,
        "work_dir": "work",
    }
    config["input_spectrum"] = {
        "path": "inputs/spectrum.dat",
        "format": "ascii_columns",
        "wavelength_unit": "Angstrom",
        "flux_unit": "erg / (cm2 s Angstrom)",
        "columns": {
            "wavelength": 0,
            "flux": 1,
        },
        "uncertainty_path": "inputs/uncertainty.dat",
    }
    config["segments"] = {
        "path": "inputs/segments.txt",
        "wavelength_reference_path": "inputs/reference.txt",
    }
    config["masking"] = {
        "exclude_regions_path": "inputs/exclude.txt",
    }
    config["outputs"] = {
        "directory": "results",
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")

    loaded = load_config(config_path)

    assert loaded.run.grid_path == config_dir / "inputs" / "grid.h5"
    assert loaded.run.work_dir == config_dir / "work"
    assert loaded.input_spectrum.path == config_dir / "inputs" / "spectrum.dat"
    assert loaded.input_spectrum.uncertainty_path == config_dir / "inputs" / "uncertainty.dat"
    assert loaded.segments.path == config_dir / "inputs" / "segments.txt"
    assert loaded.segments.wavelength_reference_path == config_dir / "inputs" / "reference.txt"
    assert loaded.masking.exclude_regions_path == config_dir / "inputs" / "exclude.txt"
    assert loaded.outputs.directory == config_dir / "results"
    assert loaded.outputs.wavelength_type is None
    assert loaded.model.normalize_parts is None


def test_load_config_preserves_symlink_path_for_relative_rebasing(tmp_path: Path) -> None:
    shared_dir = tmp_path / "shared"
    shared_dir.mkdir()
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    config_path = shared_dir / "config.json"
    symlink_path = run_dir / "config.json"

    config = _base_config()
    config["run"] = {
        "name": "toy",
        "grid_path": "inputs/grid.h5",
        "mpi": False,
        "mpi_processes": 1,
    }
    config["input_spectrum"] = {
        "path": "inputs/spectrum.dat",
        "format": "ascii_columns",
        "wavelength_unit": "Angstrom",
        "flux_unit": "erg / (cm2 s Angstrom)",
        "columns": {
            "wavelength": 0,
            "flux": 1,
        },
    }
    config["outputs"] = {
        "directory": "results",
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")
    symlink_path.symlink_to(config_path)

    loaded = load_config(symlink_path)

    assert loaded.run.grid_path == run_dir / "inputs" / "grid.h5"
    assert loaded.input_spectrum.path == run_dir / "inputs" / "spectrum.dat"
    assert loaded.outputs.directory == run_dir / "results"
    assert loaded.source_path == symlink_path


def test_load_config_rejects_unknown_keys(tmp_path: Path) -> None:
    config_path = tmp_path / "config.json"
    config = _base_config()
    config["outputs"] = {
        "directory": "out",
        "bogus": True,
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")

    with pytest.raises(ValueError, match=r"Unknown keys in 'outputs': bogus"):
        load_config(config_path)


def test_load_config_accepts_output_wavelength_type(tmp_path: Path) -> None:
    config_path = tmp_path / "config.json"
    config = _base_config()
    config["outputs"] = {
        "directory": "out",
        "wavelength_type": "air",
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")

    loaded = load_config(config_path)

    assert loaded.outputs.wavelength_type == "air"


def test_load_config_rejects_string_boolean_values(tmp_path: Path) -> None:
    config_path = tmp_path / "config.json"
    config = _base_config()
    config["run"] = {
        "name": "toy",
        "grid_path": "grid.h5",
        "mpi": "false",
        "mpi_processes": 1,
    }
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="ascii")

    with pytest.raises(ValueError, match=r"run\.mpi must be a boolean"):
        load_config(config_path)


def test_load_config_rejects_duplicate_keys(tmp_path: Path) -> None:
    config_path = tmp_path / "config.json"
    config_path.write_text(
        """
{
  "run": {
    "name": "toy",
    "grid_path": "grid.h5",
    "mpi": false,
    "mpi": true
  },
  "input_spectrum": {
    "path": "spectrum.dat",
    "format": "ascii_columns",
    "wavelength_unit": "Angstrom",
    "flux_unit": "erg / (cm2 s Angstrom)"
  },
  "model": {},
  "priors": {
    "teff": {
      "type": "uniform",
      "min": 3000.0,
      "max": 5000.0
    }
  },
  "outputs": {
    "directory": "out"
  }
}
""".strip()
        + "\n",
        encoding="ascii",
    )

    with pytest.raises(ValueError, match=r"Duplicate JSON key\(s\) detected: mpi"):
        load_config(config_path)
