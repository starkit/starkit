from __future__ import annotations

import builtins
import importlib
import sys
from types import SimpleNamespace

import pytest

from starkit.analysis.config import ModelConfig, PriorConfig


def test_pipeline_import_does_not_eagerly_import_multinest(monkeypatch: pytest.MonkeyPatch) -> None:
    sys.modules.pop("starkit.analysis.pipeline", None)
    real_import = builtins.__import__
    attempted_imports: list[str] = []

    def guarding_import(name, globals=None, locals=None, fromlist=(), level=0):
        if name == "starkit.fitkit.multinest.base":
            attempted_imports.append(name)
            raise AssertionError("MultiNest was imported while loading starkit.analysis.pipeline")
        return real_import(name, globals, locals, fromlist, level)

    monkeypatch.setattr(builtins, "__import__", guarding_import)

    module = importlib.import_module("starkit.analysis.pipeline")

    assert module.__name__ == "starkit.analysis.pipeline"
    assert attempted_imports == []


def test_configure_fit_priors_orders_by_model_parameter_order() -> None:
    pipeline = importlib.import_module("starkit.analysis.pipeline")

    class FakeFitModel:
        def __init__(self) -> None:
            self.param_names = ["teff", "mh", "logg"]
            self.fixed: dict[str, bool] = {}

    fit_model = FakeFitModel()
    prior_objects, user_names, actual_names = pipeline._configure_fit_priors(
        fit_model,
        {
            "mh": PriorConfig(type="fixed", value=-0.5),
            "teff": PriorConfig(type="uniform", min=3500.0, max=4500.0),
        },
        {
            "teff": "teff",
            "mh": "mh",
            "logg": "logg",
        },
    )

    assert user_names == ["teff", "mh"]
    assert actual_names == ["teff", "mh"]
    assert fit_model.fixed == {"teff": False, "mh": False, "logg": True}
    assert prior_objects[0](0.25) == pytest.approx(3750.0)
    assert prior_objects[1](0.25) == pytest.approx(-0.5)


def test_build_normalizer_uses_global_normalize_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    pipeline = importlib.import_module("starkit.analysis.pipeline")
    calls: list[tuple[str, object, object]] = []

    def fake_normalize(observed, npol):
        calls.append(("Normalize", observed, npol))
        return "global-normalizer"

    monkeypatch.setattr(pipeline, "Normalize", fake_normalize)

    normalizer = pipeline._build_normalizer(
        SimpleNamespace(model=ModelConfig(normalize_npol=2)),
        observed="observed-spectrum",
    )

    assert normalizer == "global-normalizer"
    assert calls == [("Normalize", "observed-spectrum", 2)]


def test_build_normalizer_uses_normalize_parts_for_segments(monkeypatch: pytest.MonkeyPatch) -> None:
    pipeline = importlib.import_module("starkit.analysis.pipeline")
    calls: list[tuple[str, object, object, object]] = []

    def fake_normalize_parts(observed, parts, npol):
        calls.append(("NormalizeParts", observed, parts, npol))
        return "parts-normalizer"

    monkeypatch.setattr(pipeline, "NormalizeParts", fake_normalize_parts)

    parts = ["part-a", "part-b"]
    normalizer = pipeline._build_normalizer(
        SimpleNamespace(model=ModelConfig(normalize_npol=3, normalize_parts="segments")),
        observed="observed-spectrum",
        normalization_parts=parts,
    )

    assert normalizer == "parts-normalizer"
    assert calls == [("NormalizeParts", "observed-spectrum", parts, 3)]


def test_build_normalizer_uses_normalize_parts_for_autodetect(monkeypatch: pytest.MonkeyPatch) -> None:
    pipeline = importlib.import_module("starkit.analysis.pipeline")
    calls: list[tuple[str, object, object, object]] = []

    def fake_normalize_parts(observed, parts, npol):
        calls.append(("NormalizeParts", observed, parts, npol))
        return "autodetect-normalizer"

    monkeypatch.setattr(pipeline, "NormalizeParts", fake_normalize_parts)

    parts = ["autodetected-part"]
    normalizer = pipeline._build_normalizer(
        SimpleNamespace(model=ModelConfig(normalize_npol=2, normalize_parts="autodetect")),
        observed="observed-spectrum",
        normalization_parts=parts,
    )

    assert normalizer == "autodetect-normalizer"
    assert calls == [("NormalizeParts", "observed-spectrum", parts, 2)]
