******************************************
Building a Grid from the PHOENIX Library
******************************************

StarKit builds grids from the high-resolution spectra of the
PHOENIX-ACES-AGSS-COND-2011 library (`Husser et al. 2013, A&A 553, A6
<https://doi.org/10.1051/0004-6361/201219058>`_). The grids it writes record
this reference in their ``bibtex`` metadata entry. Building a grid needs the
``grid_io`` extra (``pip install -e '.[grid_io]'``).

Downloading the spectra
=======================

The spectra are listed at https://phoenix.astro.physik.uni-goettingen.de/spectra/
and can be downloaded from https://ftp.gwdg.de/pub/misc/phoenix/HiResFITS/.
Download the wavelength file ``WAVE_PHOENIX-ACES-AGSS-COND-2011.fits`` and
the parts of the ``PHOENIX-ACES-AGSS-COND-2011`` directory that are needed,
keeping its layout: one directory per metallicity, such as ``Z-0.0``, and one
per metallicity and alpha enhancement, such as ``Z-0.0.Alpha=+0.20``. The
wavelength file and the ``PHOENIX-ACES-AGSS-COND-2011`` directory must be in
the same directory.

Indexing the spectra
====================

From the directory that contains the wavelength file, ``make_grid_info``
reads the parameters of every spectrum from its FITS header (``PHXTEFF``,
``PHXLOGG``, ``PHXM_H`` and ``PHXALPHA``) and writes them, with the file
names, the wavelengths and the grid metadata, to an HDF5 file::

    from starkit.gridkit.io.phoenix import make_grid_info

    make_grid_info('phoenix_grid_info.h5')

The file names are stored relative to this directory, so the next steps must
also be run from it. The index is read back with::

    import pandas as pd
    from astropy import units as u

    meta = pd.read_hdf('phoenix_grid_info.h5', 'meta')
    raw_index = pd.read_hdf('phoenix_grid_info.h5', 'index')
    wavelength = (pd.read_hdf('phoenix_grid_info.h5', 'wavelength')[0].values
                  * u.Unit(meta['wavelength_unit']))

``raw_index`` has the columns ``teff``, ``logg``, ``mh``, ``alpha`` and
``filename``, with one row per spectrum.

Selecting the spectra
=====================

Select the spectra of the grid with pandas::

    index = raw_index[raw_index.teff.between(2000, 20000)
                      & raw_index.logg.between(-1, 5)
                      & raw_index.mh.between(-1.5, 0.5)
                      & (raw_index.alpha == 0.0)]

A parameter that has a single value in the selection, such as ``alpha``
here, is not a parameter of the resulting grid.

Writing the grid
================

``PhoenixProcessGrid`` reduces each spectrum to resolving power ``R`` with
``R_sampling`` pixels per resolution element, on wavelengths evenly spaced in
their logarithm from ``wavelength_start`` to ``wavelength_stop``::

    from starkit.gridkit.io.phoenix import PhoenixProcessGrid

    pgrid = PhoenixProcessGrid(index, wavelength, meta,
                               wavelength_start=2000 * u.angstrom,
                               wavelength_stop=9000 * u.angstrom,
                               R=10000, R_sampling=4, pre_sampling=4)
    pgrid.to_hdf('phoenix_r10000.h5')

Each spectrum is first interpolated onto wavelengths evenly spaced in their
logarithm, with ``pre_sampling`` pixels per resolution element of the
PHOENIX spectra, taken as :math:`R = 30000` from 500 to 3000 Å,
:math:`5 \times 10^5` from 3000 to 25000 Å and :math:`10^5` from 25000 to
55000 Å. It is then convolved with a Gaussian that brings it to resolving
power ``R`` and interpolated onto the output wavelengths.

Set ``pre_sampling`` equal to ``R_sampling``. The width of the Gaussian is
computed with ``R_sampling`` pixels per resolution element, and with the
defaults, ``R_sampling=4`` and ``pre_sampling=2``, it is twice the intended
width. Processed with ``R=10000`` and the defaults, a line at 5000 Å with an
intrinsic FWHM of 0.01 Å has a FWHM of 1.000 Å (from a Gaussian fit), a
resolving power of 5000; with ``pre_sampling=4`` its FWHM is 0.500 Å, a
resolving power of 10000. The metadata of the grid records ``R`` in both
cases.

The fluxes are converted from the erg s\ :sup:`-1` cm\ :sup:`-2` cm\ :sup:`-1`
of the PHOENIX files to erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1`, and the
wavelengths are vacuum wavelengths in Å. Because the fluxes are surface
fluxes, the ``Distance`` operation cannot be used with these grids (see
:doc:`flux_units`). ``to_hdf`` stops with an error if the file exists,
unless ``overwrite=True`` is given.

Checking the grid
=================

``inspect_grid`` (see :doc:`gridkit_fileformat`) reports the parameters,
their sampling and the wavelength range of the new grid without reading the
fluxes, and ``load_grid`` loads it::

    from starkit.gridkit import load_grid

    grid = load_grid('phoenix_r10000.h5', wavelength_type='vacuum')

The notebook :doc:`reading_phoenix` shows the same steps.
