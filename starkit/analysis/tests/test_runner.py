from __future__ import annotations

import builtins
import json
from pathlib import Path
import subprocess

import pytest

from starkit.analysis import runner


def test_run_analysis_checks_mpi_before_importing_pipeline(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    config_path = tmp_path / "config.json"
    config_path.write_text(
        json.dumps(
            {
                "run": {
                    "name": "toy",
                    "grid_path": "/tmp/grid.h5",
                    "mpi": True,
                    "mpi_processes": 2,
                },
                "input_spectrum": {
                    "path": "/tmp/spec.dat",
                    "format": "ascii_columns",
                    "wavelength_unit": "Angstrom",
                    "flux_unit": "erg / (cm2 s Angstrom)",
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
                    "directory": str(tmp_path),
                },
            }
        )
        + "\n",
        encoding="ascii",
    )

    real_import = builtins.__import__
    attempted_imports: list[str] = []

    def guarding_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "starkit.analysis.pipeline":
            attempted_imports.append(name)
            raise AssertionError("pipeline was imported before MPI launch handling")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarding_import)

    def stop_before_pipeline(*_args, **_kwargs):
        raise RuntimeError("launch checked")

    monkeypatch.setattr(runner, "_maybe_launch_under_mpi", stop_before_pipeline)

    with pytest.raises(RuntimeError, match="launch checked"):
        runner.run_analysis(str(config_path))

    assert attempted_imports == []


def test_format_missing_file_error_mentions_config_directory(tmp_path: Path) -> None:
    config_dir = tmp_path / "cfg"
    config_dir.mkdir()
    config_path = (config_dir / "config.json").resolve()
    missing_path = str((config_dir / "inputs" / "missing.dat").resolve())

    exc = FileNotFoundError(2, "No such file or directory", missing_path)

    message = runner._format_missing_file_error(
        "Could not load input spectrum",
        exc,
        config_path,
    )

    assert "Could not load input spectrum" in message
    assert missing_path in message
    assert str(config_path) in message
    assert str(config_dir.resolve()) in message
    assert "not the current working directory" in message


def test_mpi_launcher_uses_stable_module_command(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    for key in runner._MPI_ENV_KEYS:
        monkeypatch.delenv(key, raising=False)

    monkeypatch.setattr(runner.shutil, "which", lambda name: "/usr/bin/mpiexec" if name == "mpiexec" else None)
    calls: list[list[str]] = []
    monkeypatch.setattr(subprocess, "call", lambda command: calls.append(command) or 0)

    with pytest.raises(SystemExit) as exc_info:
        runner._maybe_launch_under_mpi("/tmp/analysis.json", True, 3)

    assert exc_info.value.code == 0
    assert calls == [
        [
            "/usr/bin/mpiexec",
            "-n",
            "3",
            runner.sys.executable,
            "-m",
            "starkit.analysis",
            "run",
            "/tmp/analysis.json",
        ]
    ]
