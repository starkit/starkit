from __future__ import annotations

import json
from pathlib import Path

import pytest

from starkit.analysis import cli
from starkit.analysis import runner


def test_run_subcommand_dispatches_to_runner(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(runner, "run_analysis", lambda path: calls.append(path) or 0)

    assert cli.main(["run", "analysis.json"]) == 0
    assert calls == ["analysis.json"]


def test_validate_subcommand_can_emit_json(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    config_path = tmp_path / "invalid.json"
    config_path.write_text("{}\n", encoding="ascii")

    assert cli.main(["validate", str(config_path), "--shallow", "--json"]) == 1
    report = json.loads(capsys.readouterr().out)
    assert report["ok"] is False
    assert report["errors"]


def test_inspect_grid_subcommand_writes_requested_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from starkit.analysis import templates
    from starkit.gridkit import inspector

    inspection = object()
    monkeypatch.setattr(inspector, "inspect_grid", lambda path: inspection)
    monkeypatch.setattr(inspector, "format_inspection_summary", lambda value: "grid summary")
    monkeypatch.setattr(inspector, "inspection_to_dict", lambda value: {"path": "/tmp/grid.h5"})
    monkeypatch.setattr(
        templates,
        "build_template_config",
        lambda value, **kwargs: {"run": {"name": kwargs["run_name"]}},
    )

    summary_path = tmp_path / "summary.json"
    template_path = tmp_path / "analysis.json"
    result = cli.main(
        [
            "inspect-grid",
            "/tmp/grid.h5",
            "--summary-json-out",
            str(summary_path),
            "--template-out",
            str(template_path),
            "--name",
            "toy",
        ]
    )

    assert result == 0
    assert json.loads(summary_path.read_text(encoding="ascii")) == {"path": "/tmp/grid.h5"}
    assert json.loads(template_path.read_text(encoding="ascii")) == {"run": {"name": "toy"}}
    assert "grid summary" in capsys.readouterr().out
