"""Regression tests for the whitespace-separated spectrum readers used to build grids.

pandas 3 removed ``delim_whitespace`` and pandas 2 rejects ``header=-1``,
which these readers used before.
"""

import bz2
import importlib

import numpy as np
import pytest

pytest.importorskip("progressbar", reason="needs the grid_io extra")

# Irregular whitespace: leading spaces, tabs and repeated spaces.
ROWS = "  10000.0    1.5\n10001.0\t2.5\n   10002.0  3.5  \n"
FLUX = [1.5, 2.5, 3.5]


def test_bosz_reader_converts_bz2_text_to_npy(tmp_path):
    from starkit.gridkit.io.bosz.base import convert_bz2_memmap

    fname = tmp_path / "spectrum.asc.bz2"
    fname.write_bytes(bz2.compress(ROWS.encode()))
    convert_bz2_memmap(str(fname))
    np.testing.assert_array_equal(np.load(tmp_path / "spectrum.asc.v1.npy"), FLUX)


def test_cmfgen_reader_skips_header_row(tmp_path):
    from starkit.gridkit.io.cmfgen.base import convert_bz2_memmap

    fname = tmp_path / "spectrum.tsv"
    fname.write_text("wavelength flux\n" + ROWS)
    convert_bz2_memmap(str(fname))
    np.testing.assert_array_equal(np.load(tmp_path / "spectrum.v1.npy"), FLUX)


def test_pollux_reader_converts_headerless_text_to_npy(tmp_path):
    from starkit.gridkit.io.pollux.base import convert_bz2_memmap

    fname = tmp_path / "model.spec"
    fname.write_text(ROWS)
    convert_bz2_memmap(str(fname))
    np.testing.assert_array_equal(np.load(tmp_path / "model.spec.v1.npy"), FLUX)


def test_pollux_meta_script_does_not_build_grid_on_import():
    # The script builds the grid from a hard-coded /Volumes path; importing it
    # must not start that.
    module = importlib.import_module("starkit.gridkit.io.pollux.mk_pollux_cmfgen_meta")
    assert callable(module.main)
