***************
Getting Started
***************

.. _working-with-grid:

Loading a grid
==============

StarKit computes a model spectrum by interpolating linearly between the
spectra of a model grid stored in an HDF5 file. :doc:`io/available_grids`
lists grids that can be downloaded, and :ref:`io` describes how to build one.
The examples on this page use ``phoenix_t4000_10000_w3000_9000_r3000.h5``, a
PHOENIX grid with resolving power :math:`R = 3000` from 3000 to 9000 Å, which
is also in the ``docs/examples`` directory of the repository::

    >>> from starkit.gridkit import load_grid
    >>> grid = load_grid('phoenix_t4000_10000_w3000_9000_r3000.h5',
    ...                  wavelength_type='vacuum')
    >>> grid
    <SpectralGrid(teff=4000., logg=3., mh=0.)>

The grid is an `astropy.modeling
<https://docs.astropy.org/en/stable/modeling/>`_ model with one parameter per
grid dimension, here ``teff`` (K), ``logg`` (log g, with g in
cm s\ :sup:`-2`) and ``mh`` ([M/H], dex). The
parameters start at the values of the first spectrum in the grid file.
``wavelength_type`` is ``'air'`` or ``'vacuum'``; if it differs from the
convention of the grid, StarKit converts the wavelengths. Without it, the
convention of the grid is used and a warning is logged.

The range of each parameter is::

    >>> grid.param_names
    ('teff', 'logg', 'mh')
    >>> grid.get_grid_extent()
    [(4000.0, 10000.0), (3.0, 5.0), (-1.0, 0.0)]

Evaluating the grid
===================

Set the parameters and call the grid to obtain the spectrum::

    >>> grid.teff = 5780.
    >>> grid.logg = 4.4
    >>> grid.mh = 0.0
    >>> wavelength, flux = grid()

``wavelength`` and ``flux`` are numpy arrays in the units of the grid, given by
``grid.wavelength.unit`` and ``grid.fluxes.unit`` (Å and
erg s\ :sup:`-1` Å\ :sup:`-1` for this grid). Parameter values outside the grid
are accepted, but the flux is then NaN at every wavelength.

To plot the spectrum::

    >>> import matplotlib.pyplot as plt
    >>> plt.plot(wavelength, flux)

The notebook :doc:`/examples/interactive_spectrum` plots the spectrum with
sliders for the grid parameters.

Operations on the model spectrum
================================

Operations take a wavelength and a flux array and return modified arrays.
They are also astropy models, so their parameters can be fitted together with
the grid parameters. The operations in ``starkit.base.operations`` are listed
below; import each one from its module by name, for example
``from starkit.base.operations.stellar import RadialVelocity``, since
``RadialVelocity`` and ``Distance`` are not imported by
``from starkit.base.operations.stellar import *``.

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Operation
     - Parameters
     - Effect
   * - ``stellar.RotationalBroadening``
     - ``vrot`` (km/s), ``limb_darkening`` (fixed at 0.6)
     - Convolves with the rotational broadening profile of Gray, *The
       Observation and Analysis of Stellar Photospheres*. Create it with
       ``RotationalBroadening.from_grid(grid, vrot=...)``, which reads the
       velocity width of a grid pixel, :math:`c / (R \times R_\mathrm{sampling})`,
       from the grid. Rotation slower than one pixel width (24.98 km/s for the
       test grid) leaves the spectrum unchanged.
   * - ``stellar.DopplerShift``
     - ``vrad`` (km/s)
     - Multiplies the wavelengths by :math:`\sqrt{(1 + \beta) / (1 - \beta)}`,
       with :math:`\beta = v_\mathrm{rad} / c`.
   * - ``stellar.RadialVelocity``
     - ``vrad`` (km/s)
     - Multiplies the wavelengths by :math:`1 + \beta`.
   * - ``stellar.CCM89Extinction``
     - ``a_v`` (mag), ``r_v`` (fixed at 3.1)
     - Applies the Cardelli, Clayton & Mathis (1989) extinction law between
       0.091 and 3.33 µm; the flux outside this range is not changed.
   * - ``stellar.Distance``
     - ``distance`` (pc)
     - Divides the flux by :math:`4 \pi d^2`, which is correct only for a
       grid of luminosity density (erg s\ :sup:`-1` Å\ :sup:`-1`). Many
       grids, including the test grid above, hold surface flux, and with them
       ``Distance`` either stops with an error or returns a flux too low by
       :math:`4 \pi R^2`. Read :doc:`io/flux_units` before using it.
   * - ``spectrograph.InstrumentRConstant``
     - ``R``
     - Convolves with a Gaussian to resolving power ``R``, which must be lower
       than the resolving power of the grid; otherwise the evaluation stops
       with an ``AssertionError``. Create it with
       ``InstrumentRConstant.from_grid(grid, R=...)``.
   * - ``imager.Photometry``
     - none
     - Converts the spectrum to magnitudes in a set of filters, assuming the
       flux is in erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1` at the observer
       (see :doc:`io/flux_units`). It requires the ``photometry`` extra
       (wsynphot); see :doc:`/examples/starkit_photometry`.

For example, a Doppler shift of 20 km/s::

    >>> from starkit.base.operations.stellar import DopplerShift
    >>> doppler = DopplerShift(vrad=20.)
    >>> new_wavelength, new_flux = doppler(wavelength, flux)

Combining operations
====================

The ``|`` operator passes the output of one model to the next
(see `combining models
<https://docs.astropy.org/en/stable/modeling/compound-models.html>`_ in the
astropy documentation)::

    >>> model = grid | doppler
    >>> model
    <CompoundModel(teff_0=5780., logg_0=4.4, mh_0=0., vrad_1=20.)>

The parameter names of the combined model end in the position of the model
they belong to. Further operations are added the same way::

    >>> from starkit.base.operations.stellar import (RotationalBroadening,
    ...                                              CCM89Extinction)
    >>> rotation = RotationalBroadening.from_grid(grid, vrot=200.)
    >>> extinction = CCM89Extinction(a_v=1.0)
    >>> model = grid | rotation | doppler | extinction
    >>> model
    <CompoundModel(teff_0=5780., logg_0=4.4, mh_0=0., vrot_1=200., limb_darkening_1=0.6, vrad_2=20., a_v_3=1., r_v_3=3.1)>
    >>> model.vrot_1 = 300.
    >>> wavelength, flux = model()

Assembling a model
==================

``assemble_model`` builds the combined model from parameter names. Each
keyword adds the operation that has a parameter of that name, with that
parameter as the starting value::

    >>> from starkit import assemble_model
    >>> model = assemble_model(grid, vrot=300., vrad=20., a_v=0.3)
    >>> print(model)
    Model: CompoundModel
    Inputs: ()
    Outputs: ('wavelength', 'flux')
    Model set size: 1
    Expression: [0] | [1] | [2] | [3]
    Components:
        [0]: <SpectralGrid(teff=5780., logg=4.4, mh=0.)>

        [1]: <RotationalBroadening(vrot=300., limb_darkening=0.6)>

        [2]: <DopplerShift(vrad=20.)>

        [3]: <CCM89Extinction(a_v=0.3, r_v=3.1)>
    Parameters:
        teff_0 logg_0 mh_0 vrot_1 limb_darkening_1 vrad_2 a_v_3 r_v_3
        ------ ------ ---- ------ ---------------- ------ ----- -----
        5780.0    4.4  0.0  300.0              0.6   20.0   0.3   3.1

The keywords are ``vrot``, ``vrad``, ``a_v``, ``r_v``, ``distance`` and ``R``.
``filter_set`` and ``mag_type`` (``'vega'`` by default) add ``Photometry``.
``delta_lambda`` adds ``InstrumentDeltaLambdaConstant``, a Gaussian of
constant FWHM in wavelength, but with numpy 2.5 its evaluation stops with
``TypeError: only 0-dimensional arrays can be converted to Python scalars``.
The operations are applied in the order rotation, Doppler shift, extinction,
distance, instrumental resolution, whatever the order of the keywords.
``vrad`` selects ``DopplerShift``; to use ``RadialVelocity``, combine the
operations with ``|``. Grid parameters are not accepted as keywords and stop
``assemble_model`` with a ``ValueError``; set them on the returned model.

``assemble_model`` requires at least one stellar operation (``vrot``,
``vrad``, ``a_v``, ``r_v`` or ``distance``). Without one it stops with
``TypeError: object of type 'NoneType' has no len()``, also when ``R`` or
``spectrum`` is given. Passing ``vrad=0.`` avoids this.

Comparing a model with an observed spectrum
===========================================

An observed spectrum is given to StarKit as an ``SKSpectrum1D``, made from
wavelength, flux and uncertainty arrays::

    >>> from astropy import units as u
    >>> from starkit.fix_spectrum1d import SKSpectrum1D
    >>> spectrum = SKSpectrum1D.from_array(
    ...     wavelength=obs_wavelength * u.angstrom,
    ...     flux=obs_flux * u.Unit('erg / (s cm2 Angstrom)'),
    ...     uncertainty=obs_uncertainty)

The model is interpolated onto the observed wavelengths using their numerical
values, so the observed wavelengths must be in the wavelength unit of the grid
(Å for the grids StarKit builds), and the grid must be loaded with the
``wavelength_type`` of the observed spectrum. The interpolation is linear and
does not integrate over the observed pixels.

With ``spectrum``, ``assemble_model`` adds the interpolation onto the observed
wavelengths. With ``normalize_npol`` it also multiplies the model by the
polynomial in wavelength, of that degree, that best matches the model to the
observed flux, weighted by the uncertainties. The polynomial absorbs the
continuum shape and the flux scale, so the observed spectrum need not be flux
calibrated::

    >>> model = assemble_model(grid, spectrum=spectrum, vrad=20.,
    ...                        normalize_npol=3)
    >>> model
    <CompoundModel(teff_0=5780., logg_0=4.4, mh_0=0., vrad_1=20.)>
    >>> model_wavelength, model_flux = model()

``model_wavelength`` is the observed wavelength array, and ``model_flux`` is
the model on the observed wavelengths, scaled to the observed flux.
The interpolation and normalization have no free parameters. For spectra in
several parts, such as echelle orders or detector chips, ``normalize_parts``
takes a list of slices or index arrays into the observed spectrum, and each
part is normalized with its own polynomial; ``normalize_npol`` is then a
single degree or a list with one degree per part.

:doc:`fitting` describes how to fit the model parameters to the observed
spectrum.
