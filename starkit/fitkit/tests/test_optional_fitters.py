"""Regression tests for the optional fitting back ends."""

import importlib
import sys

import pytest
from astropy.modeling.models import Const1D

from starkit.fitkit.priors import UniformPrior


def test_jbopt_fitter_reports_missing_jbopt_when_run():
    # The module must import without a usable jbopt (the PyPI release is
    # Python 2 code) and fail only when a fit is started.
    from starkit.fitkit.jbopt import base as jbopt_base

    if jbopt_base.jbopt is not None:
        pytest.skip("a working jbopt is installed")
    fitter = jbopt_base.JBOpt(Const1D(amplitude=1.0), [UniformPrior(0.0, 2.0)])
    with pytest.raises(ImportError, match="jbopt"):
        fitter.run("unused")


def test_multinest_reports_how_to_install_pymultinest(monkeypatch):
    # A None entry in sys.modules makes "import pymultinest" fail as if the
    # package were not installed.
    monkeypatch.setitem(sys.modules, "pymultinest", None)
    for name in ("starkit.fitkit.multinest", "starkit.fitkit.multinest.base"):
        monkeypatch.delitem(sys.modules, name, raising=False)
    with pytest.raises(ModuleNotFoundError) as excinfo:
        importlib.import_module("starkit.fitkit.multinest")
    message = str(excinfo.value)
    assert "environment.yml" in message
    assert "pip install 'starkit[multinest]'" in message


def test_ipyparallel_helpers_import():
    pytest.importorskip("ipyparallel", reason="needs the parallel extra")
    try:
        module = importlib.import_module("starkit.fitkit.parallel.ipparallel")
    except SystemExit:
        # ipparallel imports the MultiNest fitter, and pymultinest calls
        # sys.exit() when it cannot load the MultiNest library.
        pytest.skip("needs the MultiNest library")
    assert callable(module.interactive)
