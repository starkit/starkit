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

Recent versions of conda solve this environment much faster than older ones.
If the solve takes more than a few seconds, update conda first::

    conda update -n base conda

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
``all``
    ``grid_io``, ``multinest`` and ``parallel``.
``test``
    pytest, for running the test suite.
``docs``
    Sphinx, its extensions and IPython, for building this documentation.

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

Synthetic photometry
====================

``Photometry`` needs wsynphot (https://github.com/starkit/wsynphot), which is
not on PyPI and therefore cannot be a StarKit extra. As of November 2026 it
does not install in a current environment:

- its setup script imports ``imp``, which Python 3.12 removed;
- with Python 3.11, the setup script fails because current setuptools no
  longer provides ``pkg_resources``;
- its code calls ``np.trapz``, which numpy 2.5 removed.

Photometry can therefore be used only where wsynphot is already installed
with numpy earlier than 2.5, until wsynphot is updated.

Running the tests
=================

From the repository::

    pytest

For an installed copy::

    pytest --pyargs starkit

Tests that need an optional dependency are skipped when it is not installed.

Building the documentation
==========================

Install the ``docs`` extra and `pandoc <https://pandoc.org/installing.html>`_,
which nbsphinx needs to convert the notebooks, then build the HTML pages from
the repository::

    pip install -e '.[docs]'
    make -C docs html

The pages are written to ``docs/_build/html``. The notebooks are included
with their saved output; they are not run during the build.

The workflow ``.github/workflows/docs.yml`` builds the documentation the same
way for every pull request, and on every push to ``master`` publishes it at
https://starkit.github.io/starkit/.
