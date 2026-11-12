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

With conda, create the environment defined in ``environment.yml`` from the
repository root::

    git clone https://github.com/starkit/starkit
    cd starkit
    conda env create -f environment.yml
    conda activate starkit

This environment uses Python 3.10 and conda-forge packages. It includes the
``grid_io``, ``multinest`` and ``parallel`` extras, the MultiNest library and
pytest, and installs StarKit from the checkout in editable mode. conda-forge
has the MultiNest library for Linux on x86-64 and for macOS, but not for Linux
on ARM or for Windows. On those platforms, remove ``multinest`` and
``pymultinest`` from ``environment.yml``, install pymultinest with pip, and
build the library as described below. To update an existing environment after
``environment.yml`` changes, run::

    conda env update -f environment.yml --prune

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

``pymultinest`` is a Python interface to the MultiNest library. The conda
environment above includes the library. With a pip installation the library
must be built separately, and its directory must be on ``LD_LIBRARY_PATH``
(Linux) or ``DYLD_LIBRARY_PATH`` (macOS).

The MultiNest library from conda-forge cuts the paths of its output files at
100 characters. Slightly longer paths give truncated file names, and longer
ones stop MultiNest with the Fortran error "File already opened in another
unit". With the default prefix, keep the ``run_dir`` given to ``MultiNest``
to 58 characters or fewer. A library built from
https://github.com/JohannesBuchner/MultiNest accepts paths of up to 1000
characters. On macOS, putting its directory on ``DYLD_LIBRARY_PATH`` makes
pymultinest in the conda environment load it instead of the conda-forge
library.

Running the tests
=================

From the repository::

    pytest

For an installed copy::

    pytest --pyargs starkit

Tests that need an optional dependency are skipped when it is not installed.
