from __future__ import annotations

from starkit.analysis.templates import build_template_config
from starkit.gridkit.inspector import GridInspection, ParameterInspection


def _inspection(*, compatible: bool = True) -> GridInspection:
    return GridInspection(
        path="/tmp/toy_grid.h5",
        parameter_names=["teff", "mh"],
        parameters=[
            ParameterInspection(
                name="teff",
                minimum=5000.0,
                maximum=5500.0,
                unique_count=2,
                unique_values_preview=[5000.0, 5500.0],
            ),
            ParameterInspection(
                name="mh",
                minimum=0.0,
                maximum=0.0,
                unique_count=1,
                unique_values_preview=[0.0],
            ),
        ],
        point_count=2,
        expected_cartesian_point_count=2,
        occupancy_fraction=1.0,
        wavelength_count=3,
        wavelength_min=1000.0,
        wavelength_max=1010.0,
        wavelength_unit="Angstrom",
        wavelength_type="vacuum",
        grid_type="log",
        wavelength_spacing="log" if compatible else "linear",
        starkit_compatible=compatible,
        flux_unit="erg / (cm2 s Angstrom)",
        resolution=20000.0,
        resolution_sampling=4.0,
        grid_shape=(2, 3),
        format_version="v1.0",
        warnings=[],
    )


def test_build_template_config_uses_inspected_grid_defaults() -> None:
    template = build_template_config(
        _inspection(),
        run_name="toy_run",
        output_directory="/tmp/output",
    )

    assert template["run"]["name"] == "toy_run"
    assert template["run"]["grid_wavelength_type"] == "vacuum"
    assert template["segments"]["resample_to_grid"] is True
    assert template["model"]["normalize_parts"] is None
    assert template["outputs"]["directory"] == "/tmp/output"
    assert template["outputs"]["wavelength_type"] is None
    assert template["outputs"]["write_corner_plot"] is True
    assert template["priors"]["teff"] == {
        "type": "uniform",
        "min": 5000.0,
        "max": 5500.0,
    }
    assert template["priors"]["mh"] == {"type": "fixed", "value": 0.0}


def test_build_template_config_disables_resampling_for_incompatible_grid() -> None:
    template = build_template_config(_inspection(compatible=False))

    assert template["segments"]["resample_to_grid"] is False
