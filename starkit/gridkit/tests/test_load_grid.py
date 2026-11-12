"""Regression tests for reading grids written in the original HDF5 layout."""

import warnings

import h5py
import numpy as np
import pandas as pd
import pytest
import tables
from astropy import units as u
from pandas.errors import PerformanceWarning

from starkit.gridkit import load_grid

WAVELENGTH = np.geomspace(20000.0, 20010.0, 5)


def _linear_flux(teff, logg):
    # Linear in both parameters, so linear interpolation reproduces it exactly.
    return 1.0 + 1e-4 * (teff - 5000.0) + 0.1 * (logg - 4.0) + 0.01 * np.arange(WAVELENGTH.size)


@pytest.fixture
def legacy_grid(tmp_path):
    """A grid whose metadata is a pickled object Series, as older grids were written."""
    path = tmp_path / "legacy_grid.h5"
    meta = pd.Series({
        "parameters": ["teff", "logg"],
        "wavelength_unit": "Angstrom",
        "wavelength_type": "vacuum",
        "flux_unit": "erg / (cm2 s Angstrom)",
        "R": 50000.0,
        "R_sampling": 4.0,
        "grid_type": "log",
        "format_version": "v1.0",
    })
    points = [(teff, logg) for teff in (5000.0, 5500.0, 6000.0) for logg in (4.0, 4.5)]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PerformanceWarning)
        meta.to_hdf(path, key="meta")
    pd.DataFrame(points, columns=["teff", "logg"]).to_hdf(path, key="index")
    pd.DataFrame({0: WAVELENGTH}).to_hdf(path, key="wavelength")
    with h5py.File(path, "a") as fh:
        fh.create_dataset("fluxes", data=np.array([_linear_flux(*p) for p in points]))
    return path


def _store_attributes_as_bytes(path):
    """Rewrite the pandas attributes the way pandas 0.15 under Python 2 stored them.

    Grids written then, including the example grid in docs/examples, store every
    string attribute as bytes and have no encoding attribute value.
    """
    with tables.open_file(path, "a") as fh:
        for node in fh.walk_nodes("/"):
            attrs = node._v_attrs
            for name in attrs._v_attrnamesuser:
                value = getattr(attrs, name)
                if name == "encoding":
                    setattr(attrs, name, None)
                elif isinstance(value, str):
                    setattr(attrs, name, np.bytes_(value.encode()))


def test_grid_with_bytes_attributes(legacy_grid):
    # pandas 3 raises TypeError on these attributes in read_hdf, which is why
    # starkit requires pandas < 3.
    _store_attributes_as_bytes(legacy_grid)
    grid = load_grid(str(legacy_grid), wavelength_type="vacuum")
    grid.teff, grid.logg = 5250.0, 4.25
    _, flux = grid()
    np.testing.assert_allclose(flux, [1.05, 1.06, 1.07, 1.08, 1.09], rtol=1e-12)


def test_legacy_grid_parameters_and_wavelengths(legacy_grid):
    grid = load_grid(str(legacy_grid), wavelength_type="vacuum")
    assert list(grid.param_names) == ["teff", "logg"]
    assert grid.wavelength.unit == u.AA
    np.testing.assert_allclose(grid.wavelength.value, WAVELENGTH)
    assert grid.fluxes.unit == u.Unit("erg / (cm2 s Angstrom)")


def test_legacy_grid_interpolates_between_nodes(legacy_grid):
    grid = load_grid(str(legacy_grid), wavelength_type="vacuum")
    grid.teff, grid.logg = 5250.0, 4.25
    _, flux = grid()
    np.testing.assert_allclose(flux, [1.05, 1.06, 1.07, 1.08, 1.09], rtol=1e-12)
