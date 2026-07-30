from __future__ import annotations

from pathlib import Path

import numpy as np

from starkit.evolutionkit.yyisochrones_galactic_center import get_isotable
from starkit.evolutionkit.yyisochrones_galactic_center import plot_teff_vs_logg
from starkit.evolutionkit.yyisochrones_galactic_center import predict_logg
from starkit.evolutionkit.yyisochrones_galactic_center import predict_logg_per_slice


def test_plot_teff_vs_logg_writes_requested_curves(tmp_path: Path) -> None:
    output_path = tmp_path / "teff_logg.png"
    fig, ax = plot_teff_vs_logg(path=output_path, show=False)

    try:
        assert output_path.exists()
        assert output_path.stat().st_size > 0
        assert [line.get_label() for line in ax.lines] == [
            "[Fe/H] = -0.5",
            "[Fe/H] = +0.0",
            "[Fe/H] = +0.5",
        ]
        assert ax.get_xlim()[0] > ax.get_xlim()[1]
        assert ax.get_ylim()[0] > ax.get_ylim()[1]
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)


def test_plot_teff_vs_logg_can_overlay_isochrone_rows(tmp_path: Path) -> None:
    output_path = tmp_path / "teff_logg_background.png"
    fig, ax = plot_teff_vs_logg(
        feh_values=(-0.6, 0.6),
        logg_max=4.5,
        include_data_points=True,
        path=output_path,
        show=False,
    )

    try:
        assert output_path.exists()
        assert output_path.stat().st_size > 0
        assert len(ax.collections) == 1
        assert [line.get_label() for line in ax.lines] == [
            "[Fe/H] = -0.6",
            "[Fe/H] = +0.6",
        ]
        assert ax.get_xlim()[0] > ax.get_xlim()[1]
        assert ax.get_ylim()[0] > ax.get_ylim()[1]
    finally:
        import matplotlib.pyplot as plt

        plt.close(fig)


def test_predict_logg_per_slice_matches_exact_giant_branch() -> None:
    iso = get_isotable()
    data = iso[(iso["logg"] < 4.5) & np.isclose(iso["FeH"], 0.6)]
    teff_all = np.power(10.0, data["logT"].to_numpy())
    peak_index = int(np.flatnonzero(np.isclose(teff_all, teff_all.max()))[-1])
    giant = data.iloc[peak_index:]
    teff = np.power(10.0, giant["logT"].to_numpy())
    truth = giant["logg"].to_numpy()
    feh = np.full(teff.shape, 0.6)

    global_pred = predict_logg(teff, feh, degree=5, logg_max=4.5)
    slice_pred = predict_logg_per_slice(teff, feh, branch="giant")

    global_rmse = np.sqrt(np.mean((global_pred - truth) ** 2))
    slice_rmse = np.sqrt(np.mean((slice_pred - truth) ** 2))

    assert slice_rmse < global_rmse
    assert slice_rmse < 1e-8


def test_predict_logg_per_slice_can_select_dwarf_branch() -> None:
    teff = np.array([5200.0])
    feh = np.array([0.6])

    giant = predict_logg_per_slice(teff, feh, branch="giant")
    dwarf = predict_logg_per_slice(teff, feh, branch="dwarf")

    assert float(dwarf[0]) > float(giant[0])
