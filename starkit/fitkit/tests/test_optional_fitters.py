"""Regression tests for the optional fitting back ends."""

import importlib

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


def test_ipyparallel_helpers_import():
    pytest.importorskip("ipyparallel", reason="needs the parallel extra")
    module = importlib.import_module("starkit.fitkit.parallel.ipparallel")
    assert callable(module.interactive)
