****************
Grid Flux Units
****************

A grid stores one of three kinds of flux:

Surface flux, :math:`F_\lambda`
    erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1`, the flux leaving the stellar
    surface. The flux observed at distance :math:`d` from a star of radius
    :math:`R` is :math:`F_\lambda (R/d)^2`.
Luminosity density, :math:`L_\lambda = 4 \pi R^2 F_\lambda`
    erg s\ :sup:`-1` Å\ :sup:`-1`. The observed flux is
    :math:`L_\lambda / (4 \pi d^2)`. The radius is then fixed for each model
    by the grid.
Continuum-normalized flux
    Dimensionless. It has no absolute flux scale.

The ``flux_unit`` entry of the grid metadata should say which, but for many
grids it does not (see `Reliability of the flux unit label`_). Check a grid
with ``check_flux_scale`` before using it for absolute fluxes.

Operations that need absolute fluxes
====================================

``Distance``
    Divides the flux by :math:`4 \pi d^2`, which is correct only for
    luminosity density. With a grid labeled as surface flux it stops with a
    ``UnitConversionError``. With a grid of surface flux labeled as
    luminosity density it runs, but the flux it returns is too low by
    :math:`4 \pi R^2`, with :math:`R` in cm: 57.0 mag for
    :math:`R = R_\odot`.
``Photometry``
    Computes magnitudes assuming that the flux is in
    erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1` at the observer. Without a
    scaling to the observer, the magnitudes of a surface-flux grid are those
    of the stellar surface. Colors do not depend on the flux scale.
``RelativePhotometryLikelihood``
    Compares magnitudes, so it needs absolute fluxes.
    ``PhotometryColorLikelihood`` compares colors and does not.

Spectral fits that normalize the model to the observed spectrum
(``normalize_npol`` in ``assemble_model``) do not depend on the flux scale.

Reliability of the flux unit label
==================================

- From December 2017 (commit ff2c6dc) to January 2025 (commit 8c5fc1f),
  StarKit labeled every grid it wrote ``erg/s/angstrom``, whatever the fluxes
  were. Since January 2025 it labels every grid ``erg/s/cm^2/angstrom``, also
  whatever the fluxes are; a grid of telluric transmission, for example, is
  labeled as surface flux.
- PHOENIX grids built with StarKit code older than 2018-10-19 (commit
  5566c62) hold luminosity density: that code multiplied by
  :math:`4 \pi R^2`, with the radius ``PHXREFF`` of each PHOENIX model.
  Grids built with later code hold surface flux.
- BOSZ grids built with code older than 2018-05-29 (commit 1e65adb) hold
  surface flux divided by :math:`\pi`; later ones hold surface flux.
- The date of a grid file does not show which code built it.
- Grids built outside StarKit can carry any label. A grid of
  continuum-normalized spectra labeled ``erg/s/angstrom`` has been seen.

The test grid ``phoenix_t4000_10000_w3000_9000_r3000.h5`` in
``docs/examples`` holds surface flux and is labeled ``erg/s/angstrom``.

Checking a grid
===============

``check_flux_scale`` integrates a few model spectra of a grid over its
wavelength range and divides the result by the integral of the surface flux
of a blackbody at the same effective temperature,
:math:`\pi B_\lambda(T_\mathrm{eff})`. It reads only the spectra it checks,
by default three spread from the lowest to the highest
:math:`T_\mathrm{eff}` of the grid::

    >>> from starkit.gridkit.flux_check import (check_flux_scale,
    ...                                         format_flux_scale_check)
    >>> check = check_flux_scale('phoenix_t4000_10000_w3000_9000_r3000.h5')
    >>> print(format_flux_scale_check(check))
    Grid: phoenix_t4000_10000_w3000_9000_r3000.h5
    Flux unit label: erg/s/angstrom (luminosity density)
    Wavelength range: 3000.0 to 8999.7 Angstrom
    Grid flux / blackbody surface flux, integrated over the wavelength range:
      teff=4000, logg=4, mh=-0.5: ratio 0.833, surface flux
      teff=6200, logg=4, mh=-0.5: ratio 0.989, surface flux
      teff=10000, logg=4, mh=-0.5: ratio 0.901, surface flux
    Fluxes: surface flux
    The fluxes do not agree with the flux unit label. Operations that need absolute fluxes, such as Distance and the magnitudes from Photometry, give wrong results with this grid.

The ratio is computed in the unit of the label, converted to
erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1` or erg s\ :sup:`-1` Å\ :sup:`-1`,
and identifies the fluxes:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Ratio
     - Fluxes
   * - 0.5 to 2
     - Surface flux.
   * - 0.2 to 0.5
     - Surface flux divided by :math:`\pi`.
   * - 0.05 to 0.2
     - Eddington flux, surface flux divided by :math:`4 \pi`.
   * - above :math:`10^{15}`
     - Luminosity density; the ratio is :math:`4 \pi R^2` in cm\ :sup:`2`,
       and the report gives the implied radius :math:`R`, which is about
       6 × 10\ :sup:`22` cm\ :sup:`2` for :math:`R_\odot`.
   * - any, with a median flux between 0.3 and 1.5 and no flux above 2
     - Continuum-normalized.

Model spectra differ from a blackbody by tens of percent, more so over a
narrow wavelength range: for the surface-flux grids checked so far the
ratios were 0.75 to 1.38, and for grids of surface flux divided by
:math:`\pi` 0.23 to 0.44. A model with non-finite or negative fluxes is reported as invalid and is
not used to classify the grid. ``check.models`` holds the ratio and the
interpretation of each model checked, ``check.interpretation`` the result for
the grid, and ``check.consistent_with_label`` whether it agrees with the
label (``None`` when that cannot be decided). ``n_models`` sets the number of
models checked, and ``teff_column`` the index column that holds
:math:`T_\mathrm{eff}`.

Absolute fluxes from a surface-flux grid
========================================

StarKit has no operation with a radius parameter yet. A surface flux can be
scaled to the observer with astropy's ``Scale`` model, applied to the flux
output while the wavelengths pass through ``Identity``::

    >>> import numpy as np
    >>> from astropy import units as u
    >>> from astropy.modeling.models import Identity, Scale
    >>> from starkit import assemble_model
    >>> from starkit.gridkit import load_grid
    >>> grid = load_grid('phoenix_t4000_10000_w3000_9000_r3000.h5',
    ...                  wavelength_type='vacuum')
    >>> grid.teff, grid.logg, grid.mh = 5780., 4.4, 0.0
    >>> radius, distance = 1.0 * u.R_sun, 10 * u.pc
    >>> dilution = Scale(float((radius / distance).to(1) ** 2))
    >>> model = assemble_model(grid, vrad=0.) | (Identity(1) & dilution)
    >>> model.param_names
    ('teff_0', 'logg_0', 'mh_0', 'vrad_1', 'factor_3')
    >>> wavelength, flux = model()
    >>> float(np.interp(5500., wavelength, flux))
    4.376138505982197e-11

For the Sun at 10 pc, this is 4.38 × 10\ :sup:`-11`
erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1` at 5500 Å, about V = 4.8 with a
V-band zero point of 3.6 × 10\ :sup:`-9`
erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1`; the absolute V magnitude of the
Sun is 4.81. The factor :math:`(R/d)^2` is a parameter of the model, but
:math:`R` and :math:`d` cannot be fitted separately from it.
