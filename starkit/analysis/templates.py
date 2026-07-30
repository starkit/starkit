"""Configuration templates built from inspected Starkit grids."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from starkit.gridkit.inspector import GridInspection


def build_template_config(
    inspection: GridInspection,
    *,
    run_name: str | None = None,
    output_directory: str | None = None,
) -> dict[str, Any]:
    """Build a starter analysis configuration from grid metadata."""
    name = run_name or Path(inspection.path).stem
    config: dict[str, Any] = {
        "run": {
            "name": name,
            "grid_path": inspection.path,
            "mpi": True,
            "mpi_processes": 1,
            "work_dir": None,
            "clean_up": True,
        },
        "input_spectrum": {
            "path": "/path/to/spectrum.dat",
            "format": "ascii_columns",
            "wavelength_type": inspection.wavelength_type,
            "columns": {
                "wavelength": 0,
                "flux": 1,
                "uncertainty": None,
            },
            "skiprows": 0,
            "wavelength_unit": inspection.wavelength_unit or "Angstrom",
            "flux_unit": inspection.flux_unit or "erg / (cm2 s Angstrom)",
            "snr": 30.0,
        },
        "segments": {
            "path": None,
            "mode": "file",
            "wavelength_reference_path": None,
            "wavelength_type": None,
            "resample_to_grid": bool(
                inspection.starkit_compatible
                and inspection.resolution is not None
                and inspection.resolution_sampling is not None
            ),
        },
        "masking": {
            "wavelength_type": None,
            "include_regions_path": None,
            "exclude_regions_path": None,
            "regions_json_path": None,
        },
        "model": {
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
        },
        "priors": {},
        "outputs": {
            "directory": output_directory or "/path/to/output",
            "wavelength_type": None,
            "write_hdf5": True,
            "write_plot": True,
            "write_corner_plot": True,
            "write_bestfit_spectrum": True,
            "write_masked_observed": True,
            "write_run_summary": True,
            "filenames": {},
        },
    }

    if inspection.wavelength_type is not None:
        config["run"]["grid_wavelength_type"] = inspection.wavelength_type

    priors: dict[str, dict[str, Any]] = {}
    for parameter in inspection.parameters:
        if np.isclose(parameter.minimum, parameter.maximum):
            priors[parameter.name] = {
                "type": "fixed",
                "value": parameter.minimum,
            }
        else:
            priors[parameter.name] = {
                "type": "uniform",
                "min": parameter.minimum,
                "max": parameter.maximum,
            }
    config["priors"] = priors
    return config
