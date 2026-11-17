"""Regression tests for reading grids written in the original HDF5 layout."""

import itertools
import json
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


def _write_legacy_grid(path, parameters, points, fluxes, extra_meta=None, json_meta=None):
    """Write a grid in the StarKit HDF5 layout.

    By default the metadata is a pickled object Series, as StarKit writes it.
    With json_meta="series" or json_meta="frame", every metadata value is
    stored as JSON text, in a Series or in a one-column DataFrame.
    """
    meta = {
        "parameters": list(parameters),
        "wavelength_unit": "Angstrom",
        "wavelength_type": "vacuum",
        "flux_unit": "erg / (cm2 s Angstrom)",
        "R": 50000.0,
        "R_sampling": 4.0,
        "grid_type": "log",
        "format_version": "v1.0",
        **(extra_meta or {}),
    }
    if json_meta is None:
        meta = pd.Series(meta)
    else:
        meta = pd.Series({key: json.dumps(value) for key, value in meta.items()})
        if json_meta == "frame":
            meta = meta.to_frame(name="value")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", PerformanceWarning)
        meta.to_hdf(path, key="meta")
    pd.DataFrame(points, columns=list(parameters)).to_hdf(path, key="index")
    pd.DataFrame({0: WAVELENGTH}).to_hdf(path, key="wavelength")
    with h5py.File(path, "a") as fh:
        fh.create_dataset("fluxes", data=np.asarray(fluxes))


@pytest.fixture
def legacy_grid(tmp_path):
    path = tmp_path / "legacy_grid.h5"
    points = [(teff, logg) for teff in (5000.0, 5500.0, 6000.0) for logg in (4.0, 4.5)]
    _write_legacy_grid(path, ["teff", "logg"], points, [_linear_flux(*p) for p in points])
    return path


def _bosz_like_points():
    """Parameter points with the BOSZ axes and a minimum logg that rises with teff.

    Qhull cannot triangulate these points in their own units (QH6271 qhull
    topology error), as with the BOSZ grids in issue #82.
    """
    teff = np.r_[np.arange(3500, 12001, 250), np.arange(12500, 20001, 500),
                 np.arange(21000, 35001, 1000)].astype(float)
    logg = np.arange(1.0, 5.01, 0.5)
    mh = np.arange(-0.5, 0.76, 0.25)
    alpha = np.arange(-0.25, 0.51, 0.25)

    def min_logg(t):
        return 1.0 if t < 8250 else 2.0 if t < 12500 else 3.0 if t < 31000 else 3.5

    return np.array([p for p in itertools.product(teff, logg, mh, alpha)
                     if p[1] >= min_logg(p[0])])


def _linear_flux_4d(teff, logg, mh, alpha):
    return (1.0 + 1e-5 * teff + 0.1 * logg + 0.2 * mh + 0.3 * alpha
            + 0.01 * np.arange(WAVELENGTH.size))


def test_grid_with_rising_minimum_logg_loads_and_interpolates(tmp_path):
    path = tmp_path / "bosz_like_grid.h5"
    points = _bosz_like_points()
    _write_legacy_grid(path, ["teff", "logg", "mh", "alpha"], points,
                       [_linear_flux_4d(*p) for p in points])
    grid = load_grid(str(path), wavelength_type="vacuum")
    grid.teff, grid.logg, grid.mh, grid.alpha = 5125.0, 3.25, 0.1, 0.05
    _, flux = grid()
    np.testing.assert_allclose(flux, _linear_flux_4d(5125.0, 3.25, 0.1, 0.05), rtol=1e-10)


def test_grid_extent_is_in_parameter_units(legacy_grid):
    grid = load_grid(str(legacy_grid), wavelength_type="vacuum")
    np.testing.assert_allclose(grid.get_grid_extent(), [(5000.0, 6000.0), (4.0, 4.5)])


def _six_point_grid(path, **kwargs):
    points = [(teff, logg) for teff in (5000.0, 5500.0, 6000.0) for logg in (4.0, 4.5)]
    _write_legacy_grid(path, ["teff", "logg"], points, [_linear_flux(*p) for p in points], **kwargs)


@pytest.mark.parametrize("layout", ["series", "frame"])
def test_grid_with_json_metadata(tmp_path, layout):
    path = tmp_path / "json_meta_grid.h5"
    _six_point_grid(path, json_meta=layout)
    grid = load_grid(str(path), wavelength_type="vacuum")
    assert list(grid.param_names) == ["teff", "logg"]
    assert grid.R == 50000.0
    assert grid.meta_grid["format_version"] == "v1.0"
    grid.teff, grid.logg = 5250.0, 4.25
    _, flux = grid()
    np.testing.assert_allclose(flux, [1.05, 1.06, 1.07, 1.08, 1.09], rtol=1e-12)


def test_text_metadata_that_parses_as_json_stays_text(tmp_path):
    # Only metadata stored entirely as JSON text is decoded. In grids written
    # by StarKit, a text value such as "1.0" or "null" stays text.
    path = tmp_path / "legacy_grid.h5"
    _six_point_grid(path, extra_meta={"grid_name": "1.0", "comment": "null"})
    grid = load_grid(str(path), wavelength_type="vacuum")
    assert grid.meta_grid["grid_name"] == "1.0"
    assert grid.meta_grid["comment"] == "null"


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
