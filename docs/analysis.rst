=================
Analysis workflow
=================

``starkit.analysis`` provides a configuration-driven workflow for preparing
observed spectra, assembling Starkit models, running MultiNest fits, and
writing results.

Command-line interface
======================

The installed command has three subcommands::

    starkit-analysis inspect-grid /path/to/grid.h5 --template-out analysis.json
    starkit-analysis validate analysis.json
    starkit-analysis run analysis.json

The equivalent module invocation is::

    python -m starkit.analysis inspect-grid /path/to/grid.h5
    python -m starkit.analysis validate analysis.json
    python -m starkit.analysis run analysis.json

``validate`` performs a deep preflight by default. It loads the grid and
spectrum and constructs the fit model. Use ``--shallow`` to check only the
configuration and referenced paths, or ``--json`` for a machine-readable
report.

Grid inspection
===============

The reusable inspection API belongs to ``starkit.gridkit.inspector``::

    from starkit.gridkit.inspector import inspect_grid

    inspection = inspect_grid("/path/to/grid.h5")

It reports interpolation parameters, grid occupancy, wavelength sampling,
resolution, metadata, and Starkit compatibility. Analysis-specific template
generation is kept separately in ``starkit.analysis.templates``::

    from starkit.analysis.templates import build_template_config

    template = build_template_config(inspection)

Configuration
=============

Configuration files contain these top-level sections:

``run``
    Run name, grid path, MPI settings, working directory, and cleanup policy.

``input_spectrum``
    Spectrum path, input format, columns, units, wavelength convention, and
    uncertainty or signal-to-noise configuration.

``segments``
    Optional wavelength segments, wavelength reference data, and grid
    resampling behavior.

``masking``
    Inclusion and exclusion files or JSON exported by the selection tool.

``model``
    Likelihood, broadening, Doppler shift, instrumental resolution,
    interpolation, cropping, and normalization.

``priors``
    Fixed, uniform, or Gaussian priors for parameters present on the assembled
    model.

``outputs``
    Output directory, wavelength convention, enabled products, and optional
    filename overrides.

Relative paths are resolved from the directory containing the configuration,
not from the process working directory. Unknown configuration keys and
duplicate JSON keys are rejected.

The packaged starter configuration is located at
``starkit/analysis/data/example_config.json``.

Model construction
==================

The configured operations are assembled in this order::

    grid
      -> rotational broadening
      -> Doppler shift
      -> instrument convolution
      -> crop to observed range
      -> interpolate
      -> normalize
      -> likelihood

Adding a model component makes its parameter available. Adding a prior makes
that parameter fitted; without a prior it remains fixed at its configured
value.

Instrument resolution must be lower than the native grid resolution because
the convolution can degrade a grid but cannot sharpen it.

MPI and MultiNest
=================

Install MPI support with::

    python -m pip install 'starkit[analysis_mpi]'

When ``run.mpi`` is true and ``run.mpi_processes`` is greater than one, the
workflow relaunches itself under ``mpiexec`` or ``mpirun``. ``mpi4py`` requires
a working MPI installation.

``pymultinest`` is a Python interface to the native MultiNest library. The
native library must be installed and visible to the runtime linker before a
fit can run.

Outputs
========

The workflow can write:

* the MultiNest result in HDF5 format;
* the best-fitting model spectrum;
* the masked observed spectrum;
* a JSON summary containing fitted and derived parameters;
* an observed-versus-model plot; and
* a posterior corner plot.

The output wavelength convention can be left in the analysis/grid convention
or converted explicitly to air or vacuum wavelengths.

YYIsochrones
============

The Galactic Center YYIsochrones implementation is available independently
from the analysis workflow::

    from starkit.evolutionkit.yyisochrones_galactic_center import predict_logg

The module preserves the Galactic Center-specific table, checksum, prediction
methods, and plotting helpers. Analysis summaries use it when a derived
surface gravity is needed.

Migration from ``starkit_analyse``
==================================

The previous scripts map to the unified interface as follows:

==============================  ==============================================
Previous command                New command
==============================  ==============================================
``starkit_analyse.py CONFIG``   ``starkit-analysis run CONFIG``
``config_lint.py CONFIG``       ``starkit-analysis validate CONFIG``
``grid_inspector.py GRID``      ``starkit-analysis inspect-grid GRID``
==============================  ==============================================

Python imports move from ``starkit_analyse`` to ``starkit.analysis``. Existing
JSON configuration structure and relative-path behavior are preserved.
