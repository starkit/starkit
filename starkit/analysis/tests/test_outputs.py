from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.ticker import ScalarFormatter
from starkit.utils.vacuumair_conversion import convert_vacuum2air

from starkit.analysis.config import (
    AnalysisConfig,
    InputSpectrumConfig,
    MaskingConfig,
    ModelConfig,
    OutputConfig,
    RunConfig,
    SegmentsConfig,
)
from starkit.analysis.outputs import _add_region_spans
from starkit.analysis.outputs import write_outputs
from starkit.evolutionkit.yyisochrones_galactic_center import predict_logg_per_slice


class FakeResult:
    def __init__(
        self,
        median: dict[str, float],
        mean: dict[str, float],
        posterior_data: pd.DataFrame | None = None,
    ) -> None:
        self.median = median
        self.mean = mean
        self.posterior_data = posterior_data
        self.parameter_names = [] if posterior_data is None else [
            column
            for column in posterior_data.columns
            if column not in {"weights", "loglikelihood", "x"}
        ]

    def to_hdf(self, path: str) -> None:
        Path(path).write_text("fake\n", encoding="ascii")


def _make_config(output_dir: Path) -> AnalysisConfig:
    return AnalysisConfig(
        run=RunConfig(name="toy", grid_path=Path("/tmp/grid.h5"), mpi=False, mpi_processes=1),
        input_spectrum=InputSpectrumConfig(
            path=Path("/tmp/spec.dat"),
            format="ascii_columns",
            wavelength_unit="Angstrom",
            flux_unit="erg / (cm2 s Angstrom)",
        ),
        segments=SegmentsConfig(),
        masking=MaskingConfig(),
        model=ModelConfig(),
        priors={"teff": SimpleNamespace(type="uniform"), "mh": SimpleNamespace(type="uniform")},
        outputs=OutputConfig(
            directory=output_dir,
            write_hdf5=False,
            write_plot=False,
            write_bestfit_spectrum=False,
            write_masked_observed=False,
            write_run_summary=True,
        ),
    )


def test_write_outputs_derives_logg_when_missing_from_grid(tmp_path: Path) -> None:
    config = _make_config(tmp_path)
    result = FakeResult(
        median={"teff": 3500.0, "mh": -0.2},
        mean={"teff": 3550.0, "mh": -0.1},
    )
    observed = SimpleNamespace(
        wavelength=np.array([10000.0, 10001.0]) * u.angstrom,
        flux=np.array([1.0, 1.1]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=None,
    )
    built_model = SimpleNamespace(parameter_name_map={"teff": "teff", "mh": "mh"})

    write_outputs(
        config=config,
        result=result,
        observed=observed,
        best_fit_wavelength=np.array([10000.0, 10001.0]),
        best_fit_flux=np.array([1.0, 1.1]),
        built_model=built_model,
    )

    summary_path = tmp_path / "toy_summary.json"
    assert summary_path.read_text(encoding="ascii").endswith("\n")
    summary = json.loads(summary_path.read_text(encoding="ascii"))
    assert "logg" in summary["median"]
    assert "logg" in summary["mean"]
    assert summary["median"]["logg"] == pytest.approx(
        float(np.asarray(predict_logg_per_slice(3500.0, -0.2, branch="giant")).reshape(-1)[0])
    )
    assert summary["mean"]["logg"] == pytest.approx(
        float(np.asarray(predict_logg_per_slice(3550.0, -0.1, branch="giant")).reshape(-1)[0])
    )
    assert summary["derived"]["logg_source"] == "yyisochrones_galactic_center"
    assert summary["derived"]["logg_from_grid"] is False


def test_write_outputs_keeps_existing_logg(tmp_path: Path) -> None:
    config = _make_config(tmp_path)
    result = FakeResult(
        median={"teff": 3500.0, "mh": -0.2, "logg": 1.5},
        mean={"teff": 3550.0, "mh": -0.1, "logg": 1.6},
    )
    observed = SimpleNamespace(
        wavelength=np.array([10000.0, 10001.0]) * u.angstrom,
        flux=np.array([1.0, 1.1]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=None,
    )

    write_outputs(
        config=config,
        result=result,
        observed=observed,
        best_fit_wavelength=np.array([10000.0, 10001.0]),
        best_fit_flux=np.array([1.0, 1.1]),
        built_model=SimpleNamespace(parameter_name_map={"teff": "teff", "mh": "mh", "logg": "logg"}),
    )

    summary = json.loads((tmp_path / "toy_summary.json").read_text(encoding="ascii"))
    assert summary["median"]["logg"] == 1.5
    assert "derived" not in summary


def test_write_outputs_orders_summary_parameters_by_grid_then_non_grid(tmp_path: Path) -> None:
    config = _make_config(tmp_path)
    result = FakeResult(
        median={"vrad_2": 15.0, "mh_0": -0.2, "teff_0": 3500.0, "vrot_1": 8.0},
        mean={"vrad_2": 14.0, "mh_0": -0.1, "teff_0": 3550.0, "vrot_1": 7.5},
    )
    observed = SimpleNamespace(
        wavelength=np.array([10000.0, 10001.0]) * u.angstrom,
        flux=np.array([1.0, 1.1]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=None,
    )
    built_model = SimpleNamespace(
        grid=SimpleNamespace(param_names=("teff", "mh")),
        parameter_name_map={"teff": "teff_0", "mh": "mh_0", "vrot": "vrot_1", "vrad": "vrad_2"},
    )

    write_outputs(
        config=config,
        result=result,
        observed=observed,
        best_fit_wavelength=np.array([10000.0, 10001.0]),
        best_fit_flux=np.array([1.0, 1.1]),
        built_model=built_model,
    )

    summary_text = (tmp_path / "toy_summary.json").read_text(encoding="ascii")
    summary = json.loads(summary_text)
    assert list(summary["median"]) == ["teff_0", "mh_0", "vrad_2", "vrot_1", "logg"]
    assert list(summary["mean"]) == ["teff_0", "mh_0", "vrad_2", "vrot_1", "logg"]


def test_write_outputs_writes_corner_plot_when_posterior_available(
    tmp_path: Path,
    monkeypatch,
) -> None:
    config = _make_config(tmp_path)
    config = AnalysisConfig(
        run=config.run,
        input_spectrum=config.input_spectrum,
        segments=config.segments,
        masking=config.masking,
        model=config.model,
        priors=config.priors,
        outputs=OutputConfig(
            directory=tmp_path,
            write_hdf5=False,
            write_plot=False,
            write_corner_plot=True,
            write_bestfit_spectrum=False,
            write_masked_observed=False,
            write_run_summary=False,
        ),
    )
    posterior_data = pd.DataFrame(
        {
            "weights": [0.2, 0.3, 0.5],
            "loglikelihood": [3.0, 2.0, 1.0],
            "teff_0": [3500.0, 3550.0, 3600.0],
            "mh_0": [-0.2, -0.1, 0.0],
            "vrad_1": [-10.0, 0.0, 12.0],
            "cfe_0": [0.0, 0.0, 0.0],
        }
    )
    result = FakeResult(
        median={"teff_0": 3550.0, "mh_0": -0.1},
        mean={"teff_0": 3550.0, "mh_0": -0.1},
        posterior_data=posterior_data,
    )
    observed = SimpleNamespace(
        wavelength=np.array([10000.0, 10001.0]) * u.angstrom,
        flux=np.array([1.0, 1.1]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=None,
    )
    captured: dict[str, object] = {}

    class FakeCornerModule:
        def corner(self, samples, **kwargs):
            captured["samples"] = samples
            captured["kwargs"] = kwargs
            figure = kwargs["fig"]
            axes = np.asarray(figure.subplots(samples.shape[1], samples.shape[1]), dtype=object)
            axes[0, 0].set_title(r"Teff = ${3.55e+03}_{-5.00e+01}^{+5.00e+01}$")
            axes[1, 1].set_title(r"[M/H] = ${-1.00e-01}_{-5.00e-02}^{+5.00e-02}$")
            for axis in axes.reshape(-1):
                axis.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
                axis.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
            captured["figure"] = figure
            return figure

    import starkit.analysis.outputs as outputs_module

    monkeypatch.setattr(outputs_module, "_load_corner_module", lambda: FakeCornerModule())

    write_outputs(
        config=config,
        result=result,
        observed=observed,
        best_fit_wavelength=np.array([10000.0, 10001.0]),
        best_fit_flux=np.array([1.0, 1.1]),
        built_model=SimpleNamespace(
            grid=SimpleNamespace(param_names=("teff", "mh")),
            parameter_name_map={"teff": "teff_0", "mh": "mh_0", "cfe": "cfe_0", "vrad": "vrad_1"},
        ),
    )

    assert (tmp_path / "toy_corner.pdf").exists()
    np.testing.assert_allclose(captured["samples"], posterior_data[["teff_0", "mh_0"]].to_numpy())
    assert captured["kwargs"]["labels"] == ["Teff", "[M/H]"]
    np.testing.assert_allclose(captured["kwargs"]["weights"], posterior_data["weights"].to_numpy())
    assert captured["kwargs"]["plot_contours"] is False
    assert captured["kwargs"]["fill_contours"] is False
    assert captured["kwargs"]["quiet"] is True
    figure = captured["figure"]
    teff_axis = figure.axes[0]
    assert teff_axis.get_title().startswith("Teff = ${")
    assert "e+" not in teff_axis.get_title().lower()
    assert "e-" not in teff_axis.get_title().lower()
    assert teff_axis.xaxis.get_major_formatter()(3550.4, 0) == "3550"
    assert teff_axis.yaxis.get_major_formatter()(3550.6, 0) == "3551"
    assert figure.axes[3].get_title() == r"[M/H] = ${-1.00e-01}_{-5.00e-02}^{+5.00e-02}$"


def test_write_outputs_plot_marks_include_and_exclude_regions(
    tmp_path: Path,
    monkeypatch,
) -> None:
    config = _make_config(tmp_path)
    config = AnalysisConfig(
        run=config.run,
        input_spectrum=config.input_spectrum,
        segments=config.segments,
        masking=config.masking,
        model=config.model,
        priors=config.priors,
        outputs=OutputConfig(
            directory=tmp_path,
            write_hdf5=False,
            write_plot=True,
            write_corner_plot=False,
            write_bestfit_spectrum=False,
            write_masked_observed=False,
            write_run_summary=False,
        ),
    )
    observed = SimpleNamespace(
        wavelength=np.array([10000.0, 10000.5, 10001.0]) * u.angstrom,
        flux=np.array([1.0, 1.1, 0.9]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=None,
    )
    figure, axis = plt.subplots()
    captured: dict[str, object] = {}

    import starkit.analysis.outputs as outputs_module

    def fake_subplots(*args, **kwargs):
        captured["figure"] = figure
        captured["axis"] = axis
        return figure, axis

    monkeypatch.setattr(outputs_module.plt, "subplots", fake_subplots)
    monkeypatch.setattr(outputs_module.plt, "close", lambda _: None)

    write_outputs(
        config=config,
        result=FakeResult(median={"teff": 3500.0}, mean={"teff": 3510.0}),
        observed=observed,
        best_fit_wavelength=np.array([10000.0, 10000.5, 10001.0]),
        best_fit_flux=np.array([1.0, 1.05, 0.95]),
        built_model=SimpleNamespace(
            grid=SimpleNamespace(wavelength_type="vacuum"),
            parameter_name_map={"teff": "teff"},
        ),
        include_regions=np.array([[9999.8, 10000.4]]),
        exclude_regions=np.array([[10000.6, 10001.2]]),
    )

    assert (tmp_path / "toy_fit.pdf").exists()
    assert axis.get_xlabel() == "Wavelength (Angstrom, vacuum)"
    patches = axis.patches
    assert len(patches) == 2
    assert patches[0].get_label() == "included region"
    assert patches[1].get_label() == "excluded region"
    assert patches[0].get_alpha() == pytest.approx(0.12)
    assert patches[1].get_alpha() == pytest.approx(0.22)
    plt.close(figure)


def test_write_outputs_converts_spectrum_files_to_configured_wavelength_type(tmp_path: Path) -> None:
    config = _make_config(tmp_path)
    config = AnalysisConfig(
        run=config.run,
        input_spectrum=config.input_spectrum,
        segments=config.segments,
        masking=config.masking,
        model=config.model,
        priors=config.priors,
        outputs=OutputConfig(
            directory=tmp_path,
            wavelength_type="air",
            write_hdf5=False,
            write_plot=False,
            write_corner_plot=False,
            write_bestfit_spectrum=True,
            write_masked_observed=True,
            write_run_summary=False,
        ),
    )
    vacuum_wavelength = np.array([10000.0, 10001.0])
    expected_air_wavelength = convert_vacuum2air(vacuum_wavelength * u.angstrom).to_value(u.angstrom)
    observed = SimpleNamespace(
        wavelength=vacuum_wavelength * u.angstrom,
        flux=np.array([1.0, 1.1]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=None,
    )

    write_outputs(
        config=config,
        result=FakeResult(median={"teff": 3500.0}, mean={"teff": 3510.0}),
        observed=observed,
        best_fit_wavelength=vacuum_wavelength,
        best_fit_flux=np.array([1.0, 1.05]),
        built_model=SimpleNamespace(
            grid=SimpleNamespace(wavelength_type="vacuum"),
            parameter_name_map={"teff": "teff"},
        ),
    )

    best_fit_data = np.loadtxt(tmp_path / "toy_bestfit_spectrum.txt")
    masked_observed_data = np.loadtxt(tmp_path / "toy_masked_observed.txt")
    np.testing.assert_allclose(best_fit_data[:, 0], expected_air_wavelength)
    np.testing.assert_allclose(masked_observed_data[:, 0], expected_air_wavelength)


def test_write_outputs_plot_labels_configured_wavelength_type_and_converts_axis(
    tmp_path: Path,
    monkeypatch,
) -> None:
    config = _make_config(tmp_path)
    config = AnalysisConfig(
        run=config.run,
        input_spectrum=config.input_spectrum,
        segments=config.segments,
        masking=config.masking,
        model=config.model,
        priors=config.priors,
        outputs=OutputConfig(
            directory=tmp_path,
            wavelength_type="air",
            write_hdf5=False,
            write_plot=True,
            write_corner_plot=False,
            write_bestfit_spectrum=False,
            write_masked_observed=False,
            write_run_summary=False,
        ),
    )
    vacuum_wavelength = np.array([10000.0, 10001.0])
    expected_air_wavelength = convert_vacuum2air(vacuum_wavelength * u.angstrom).to_value(u.angstrom)
    observed = SimpleNamespace(
        wavelength=vacuum_wavelength * u.angstrom,
        flux=np.array([1.0, 1.1]) * (u.erg / (u.cm**2 * u.s * u.angstrom)),
        uncertainty=None,
    )
    figure, axis = plt.subplots()
    captured: dict[str, object] = {}

    import starkit.analysis.outputs as outputs_module

    def fake_subplots(*args, **kwargs):
        captured["figure"] = figure
        captured["axis"] = axis
        return figure, axis

    monkeypatch.setattr(outputs_module.plt, "subplots", fake_subplots)
    monkeypatch.setattr(outputs_module.plt, "close", lambda _: None)

    write_outputs(
        config=config,
        result=FakeResult(median={"teff": 3500.0}, mean={"teff": 3510.0}),
        observed=observed,
        best_fit_wavelength=vacuum_wavelength,
        best_fit_flux=np.array([1.0, 1.05]),
        built_model=SimpleNamespace(
            grid=SimpleNamespace(wavelength_type="vacuum"),
            parameter_name_map={"teff": "teff"},
        ),
        include_regions=np.array([[vacuum_wavelength[0], vacuum_wavelength[1]]]),
    )

    assert axis.get_xlabel() == "Wavelength (Angstrom, air)"
    np.testing.assert_allclose(axis.lines[0].get_xdata(), expected_air_wavelength)
    np.testing.assert_allclose(axis.lines[1].get_xdata(), expected_air_wavelength)
    plt.close(figure)


def test_add_region_spans_draws_zero_width_regions_as_lines() -> None:
    figure, axis = plt.subplots()

    _add_region_spans(
        axis,
        np.array([[10000.5, 10000.5]]),
        10000.0,
        10001.0,
        color="#c92a2a",
        alpha=0.22,
        label="excluded region",
    )

    assert len(axis.patches) == 0
    assert len(axis.lines) == 1
    assert axis.lines[0].get_xdata()[0] == pytest.approx(10000.5)
    assert axis.lines[0].get_label() == "excluded region"
    plt.close(figure)
