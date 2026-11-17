*******************
GridKit File Format
*******************

A StarKit grid is an HDF5 file with four parts:

``index``
    A pandas DataFrame, written with ``to_hdf``, with one row per model and a
    column for each grid parameter.
``fluxes``
    A 2D array written with h5py, with one row per model in the order of
    ``index`` and one column per wavelength.
``wavelength``
    A pandas DataFrame whose first column holds the wavelength of each flux
    column.
``meta``
    The grid metadata, a pandas Series written with ``to_hdf``.

``load_grid`` reads these metadata entries:

``parameters``
    List of the ``index`` columns to interpolate in.
``grid_type``
    Only ``'log'`` is supported: wavelengths evenly spaced in the logarithm
    of the wavelength.
``R``
    Resolving power of the grid spectra.
``R_sampling``
    Number of pixels per resolution element.
``wavelength_unit``
    Astropy unit of ``wavelength``.
``wavelength_type``
    ``'air'`` or ``'vacuum'``.
``flux_unit``
    Astropy unit of ``fluxes``. For many grids this label is wrong; see
    :doc:`flux_units`.
``format_version``
    Grid format version, currently ``'v1.0'``.

The grids that StarKit builds also store ``bibtex`` and ``uuid``.

The metadata can be stored in two ways. StarKit writes it as a Series of
Python objects, in which ``parameters`` is a list and ``R`` is a number.
Alternatively, every entry can be stored as JSON text, in a Series or in a
one-column DataFrame: ``parameters`` as ``'["teff", "logg"]'`` and
``wavelength_type`` as ``'"vacuum"'``, for example. ``load_grid`` decodes the
metadata only when every entry is text, so text entries in the first layout,
such as a ``uuid``, are read unchanged.

Inspecting a grid
=================

``inspect_grid`` reports a grid's parameters and their sampling, the fraction
of the full parameter grid that has models, the wavelength range and sampling,
and the metadata, without reading the fluxes::

    from starkit.gridkit.inspector import inspect_grid, format_inspection_summary

    print(format_inspection_summary(inspect_grid("grid.h5")))

To check what the fluxes are, use ``check_flux_scale``, described in
:doc:`flux_units`.
