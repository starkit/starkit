************
Installation
************

StarKit requires Python 3.10 or later. Install it with pip::

    pip install git+https://github.com/starkit/starkit

For a development version, clone the repository and install it in editable
mode::

    git clone https://github.com/starkit/starkit
    cd starkit
    pip install -e '.[all,test]'

With conda, create an environment first and install StarKit into it with
pip::

    conda create -n starkit -c conda-forge python=3.12
    conda activate starkit
    pip install -e '.[all,test]'

Optional dependencies
=====================

Optional features are installed as extras, for example
``pip install -e '.[multinest]'``.

``grid_io``
    Building grids from model spectra (progressbar2, sqlalchemy).
``multinest``
    MultiNest fitting and corner plots (pymultinest, corner).
``parallel``
    Distributed fitting with ipyparallel (ipyparallel, psutil).
``photometry``
    Synthetic photometry (wsynphot). wsynphot is not on PyPI and may need to
    be installed separately.
``all``
    ``grid_io``, ``multinest`` and ``parallel``.
``test``
    pytest, for running the test suite.
``docs``
    Sphinx and its extensions, for building this documentation.

``pymultinest`` is a Python interface to the MultiNest library, which must be
built separately. Its directory must be on ``LD_LIBRARY_PATH`` (Linux) or
``DYLD_LIBRARY_PATH`` (macOS).

Running the tests
=================

From the repository::

    pytest

For an installed copy::

    pytest --pyargs starkit

Tests that need an optional dependency are skipped when it is not installed.
