from __future__ import annotations

import logging
import os
from pathlib import Path
import shutil
import subprocess
import sys

LOGGER = logging.getLogger(__name__)

_MPI_ENV_KEYS = (
    "OMPI_COMM_WORLD_SIZE",
    "PMI_SIZE",
    "PMIX_RANK",
    "MPI_LOCALNRANKS",
)


def run_analysis(config_path: str) -> int:
    from .config import load_config
    from .mpi import build_runtime_context, configure_logging

    config = load_config(config_path)
    resolved_config_path = str(config.source_path or config_path)
    _maybe_launch_under_mpi(
        resolved_config_path,
        config.run.mpi,
        config.run.mpi_processes,
    )
    runtime = build_runtime_context(config.run.mpi, config.run.work_dir)
    configure_logging(runtime.rank)

    from .pipeline import (
        build_model,
        evaluate_best_fit,
        get_grid_sampling,
        load_analysis_grid,
        run_fit,
    )
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
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            _format_missing_file_error("Could not load input spectrum", exc, config.source_path)
        ) from exc

    try:
        grid = load_analysis_grid(config)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            _format_missing_file_error("Could not load analysis grid", exc, config.source_path)
        ) from exc
    grid_wavelength_type = getattr(grid, "wavelength_type", None)
    analysis_wavelength_type = grid_wavelength_type or config.input_spectrum.wavelength_type
    try:
        analysis_arrays = prepare_analysis_arrays(
            raw_arrays,
            config.input_spectrum,
            config.segments,
            *get_grid_sampling(grid),
            grid_reference_wavelength=grid.wavelength.value,
            grid_wavelength_type=grid_wavelength_type,
        )
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            _format_missing_file_error(
                "Could not prepare analysis wavelength arrays",
                exc,
                config.source_path,
            )
        ) from exc

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
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            _format_missing_file_error("Could not apply masking", exc, config.source_path)
        ) from exc
    observed = build_observed_spectrum(masked_arrays, config.input_spectrum)
    normalization_parts = None
    if config.model.normalize_parts in {"segments", "autodetect"} and config.model.normalize_npol is not None:
        segment_bounds = masked_arrays.segment_bounds
        if config.model.normalize_parts == "autodetect":
            segment_bounds = autodetect_normalization_segment_bounds(analysis_arrays.wavelength)
        normalization_parts = build_segment_normalization_parts(
            observed,
            segment_bounds,
            config.model.normalize_npol,
        )
    built_model = build_model(config, grid, observed, normalization_parts=normalization_parts)
    result = run_fit(built_model, runtime, config.run.clean_up)

    if runtime.enabled:
        from mpi4py import MPI

        MPI.COMM_WORLD.Barrier()

    if runtime.is_root:
        from .outputs import write_outputs

        best_fit_wavelength, best_fit_flux = evaluate_best_fit(built_model.model, result)
        write_outputs(
            config=config,
            result=result,
            observed=observed,
            best_fit_wavelength=best_fit_wavelength,
            best_fit_flux=best_fit_flux,
            built_model=built_model,
            include_regions=masked_arrays.include_regions,
            exclude_regions=masked_arrays.exclude_regions,
        )

    return 0


def _maybe_launch_under_mpi(config_path: str, mpi_enabled: bool, mpi_processes: int) -> None:
    if not mpi_enabled or mpi_processes <= 1:
        return
    if any(key in os.environ for key in _MPI_ENV_KEYS):
        return

    mpiexec = shutil.which("mpiexec") or shutil.which("mpirun")
    if mpiexec is None:
        raise RuntimeError(
            "MPI is enabled and run.mpi_processes > 1, but neither 'mpiexec' nor 'mpirun' was found"
        )

    command = [
        mpiexec,
        "-n",
        str(mpi_processes),
        sys.executable,
        "-m",
        "starkit.analysis",
        "run",
        config_path,
    ]
    raise SystemExit(subprocess.call(command))


def _format_missing_file_error(context: str, exc: FileNotFoundError, config_path: Path | None) -> str:
    from .config import format_path_resolution_hint

    missing_path = exc.filename
    if missing_path is None and exc.args:
        missing_path = str(exc.args[0])

    if missing_path is None:
        return f"{context}: referenced file was not found. {format_path_resolution_hint(config_path)}"

    return (
        f"{context}: referenced file was not found: {missing_path}. "
        f"{format_path_resolution_hint(config_path)}"
    )
