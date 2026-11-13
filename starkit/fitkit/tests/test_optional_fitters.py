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


def _is_multinest_module(name):
    return name.split(".")[0] == "pymultinest" or name.startswith("starkit.fitkit.multinest")


@pytest.fixture
def fresh_multinest_import(monkeypatch):
    """Make the next import of the MultiNest fitter import pymultinest again."""
    for name in [name for name in sys.modules if _is_multinest_module(name)]:
        monkeypatch.delitem(sys.modules, name)
    yield monkeypatch
    # Drop what the test imported; monkeypatch then restores the originals.
    for name in [name for name in sys.modules if _is_multinest_module(name)]:
        del sys.modules[name]


def _write_fake_pymultinest(directory, files):
    package = directory / "pymultinest"
    package.mkdir()
    for name, code in files.items():
        (package / name).write_text(code)


def test_multinest_reports_how_to_install_pymultinest(fresh_multinest_import):
    # A None entry in sys.modules makes "import pymultinest" fail as if the
    # package were not installed.
    fresh_multinest_import.setitem(sys.modules, "pymultinest", None)
    with pytest.raises(ModuleNotFoundError) as excinfo:
        importlib.import_module("starkit.fitkit.multinest")
    message = str(excinfo.value)
    assert "environment.yml" in message
    assert "pip install 'starkit[multinest]'" in message


def test_multinest_reports_missing_library_when_pymultinest_exits(
        tmp_path, fresh_multinest_import):
    # On Linux, pymultinest calls sys.exit(1) when it cannot load the library.
    _write_fake_pymultinest(tmp_path, {"__init__.py": "import sys\nsys.exit(1)\n"})
    fresh_multinest_import.syspath_prepend(tmp_path)
    with pytest.raises(ImportError, match="could not load the MultiNest library") as excinfo:
        importlib.import_module("starkit.fitkit.multinest")
    assert "conda install -c conda-forge multinest" in str(excinfo.value)


def test_multinest_reports_missing_library_when_pymultinest_loads_none(
        tmp_path, fresh_multinest_import):
    # On macOS, pymultinest passes find_library's None to cdll.LoadLibrary,
    # which opens the Python process instead of the library. The import
    # succeeds and the fit fails later with "symbol not found".
    _write_fake_pymultinest(tmp_path, {
        "__init__.py": "from .run import run\n",
        "run.py": (
            "from ctypes import cdll\n"
            "lib = cdll.LoadLibrary(None)\n"
            "def run(*args, **kwargs):\n"
            "    lib.run(*args)\n"
        ),
    })
    fresh_multinest_import.syspath_prepend(tmp_path)
    with pytest.raises(ImportError, match="could not load the MultiNest library") as excinfo:
        importlib.import_module("starkit.fitkit.multinest")
    assert "conda install -c conda-forge multinest" in str(excinfo.value)


def test_ipyparallel_helpers_import():
    pytest.importorskip("ipyparallel", reason="needs the parallel extra")
    try:
        importlib.import_module("starkit.fitkit.multinest")
    except ImportError:
        pytest.skip("ipparallel imports the MultiNest fitter, which needs "
                    "pymultinest and the MultiNest library")
    module = importlib.import_module("starkit.fitkit.parallel.ipparallel")
    assert callable(module.interactive)
