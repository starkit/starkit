*******
Fitting
*******

StarKit fits spectra with MultiNest, through the ``MultiNest`` class in
``starkit.fitkit.multinest``. MultiNest fitting needs the pymultinest package
and the compiled MultiNest library, and a default installation of StarKit
includes neither. The conda environment in ``environment.yml`` includes both.
With pip, install the ``multinest`` extra and build the library separately.
:doc:`installation` describes both methods.

A fit has three parts: a likelihood model, which compares the model spectrum
with the observed one; a prior for each free parameter; and the MultiNest run,
which samples the posterior distribution. The example below fits the test
grid of :doc:`getting_started` to a synthetic spectrum.

The observed spectrum
=====================

The observed spectrum is an ``SKSpectrum1D``, with wavelengths in the
wavelength unit of the grid (see :doc:`getting_started`). For this example,
a spectrum with :math:`T_\mathrm{eff} = 5780` K, :math:`\log g = 4.4`,
[M/H] = -0.3 and :math:`v_\mathrm{rad} = 25` km/s is computed from the grid
between 5000 and 5500 Å, normalized to its median, and given Gaussian noise
with a signal-to-noise ratio of 50 per pixel (572 pixels)::

    >>> import numpy as np
    >>> from astropy import units as u
    >>> from starkit import assemble_model
    >>> from starkit.gridkit import load_grid
    >>> from starkit.fix_spectrum1d import SKSpectrum1D
    >>> grid = load_grid('phoenix_t4000_10000_w3000_9000_r3000.h5',
    ...                  wavelength_type='vacuum')
    >>> grid.teff, grid.logg, grid.mh = 5780., 4.4, -0.3
    >>> wavelength, flux = assemble_model(grid, vrad=25.)()
    >>> in_range = (wavelength > 5000) & (wavelength < 5500)
    >>> obs_wavelength = wavelength[in_range][::2]
    >>> obs_flux = np.interp(obs_wavelength, wavelength, flux)
    >>> obs_flux /= np.median(obs_flux)
    >>> obs_uncertainty = obs_flux / 50.
    >>> obs_flux += np.random.default_rng(0).normal(0., obs_uncertainty)
    >>> spectrum = SKSpectrum1D.from_array(
    ...     wavelength=obs_wavelength * u.angstrom,
    ...     flux=obs_flux * u.Unit('erg / (s cm2 Angstrom)'),
    ...     uncertainty=obs_uncertainty)

The likelihood
==============

The likelihood is the model followed by a likelihood operation, which returns
the natural logarithm of the likelihood::

    >>> from starkit.fitkit.likelihoods import SpectralChi2Likelihood
    >>> grid.teff, grid.logg, grid.mh = 6000., 4.0, -0.5
    >>> model = assemble_model(grid, spectrum=spectrum, vrad=0.,
    ...                        normalize_npol=3)
    >>> likelihood = model | SpectralChi2Likelihood(spectrum)
    >>> likelihood
    <CompoundModel(teff_0=6000., logg_0=4., mh_0=-0.5, vrad_1=0.)>
    >>> likelihood()
    -642.5948937184828

The normalization emits ``UserWarning: The fit may be poorly conditioned``
whenever its polynomial fit has full rank, which is the normal case, so the
warning can be ignored here.

The likelihood operations in ``starkit.fitkit.likelihoods`` for spectra are:

``SpectralChi2Likelihood(observed)``
    :math:`\ln L = -\chi^2 / 2`, with the observed uncertainties.
``SpectralChi2LikelihoodAddErr(observed, add_err=0.0)``
    Gaussian likelihood, including its normalization, with the fitted
    parameter ``add_err`` added in quadrature to the observed uncertainties.
``SpectralScaledChi2Likelihood(observed, lnf=0.0)``
    Gaussian likelihood whose variance is the observed variance plus
    :math:`f^2` times the squared observed flux, with :math:`\ln f` fitted as
    the parameter ``lnf``.
``SpectralL1Likelihood(observed)``
    :math:`\ln L = -\sum |F_\mathrm{obs} - F_\mathrm{model}| / \sigma`, for
    noise with a Laplace distribution or to reduce the weight of outliers.

If an observed spectrum has no uncertainties, all of them are taken to be 1.
A likelihood that evaluates to NaN, for example because a parameter is
outside the grid, is replaced by :math:`-10^{300}`.

Priors
======

MultiNest needs a prior for each free parameter, in the order of
``likelihood.param_names``. ``grid.get_grid_uniform_priors()`` returns uniform
priors over the range of each grid parameter::

    >>> from starkit.fitkit.priors import UniformPrior
    >>> priors = grid.get_grid_uniform_priors() + [UniformPrior(-100., 100.)]

``starkit.fitkit.priors`` also has ``GaussianPrior(m, sigma)`` and
``PoissonPrior(m)``.

Parameters marked as fixed are held at their current value and need no prior.
For example, after ``likelihood.vrad_1 = 25.`` and
``likelihood.vrad_1.fixed = True``, the radial velocity stays at 25 km/s and
the list of priors has three entries, for ``teff_0``, ``logg_0`` and
``mh_0``. The example below leaves all four parameters free.

Running the fit
===============

::

    >>> from starkit.fitkit.multinest import MultiNest
    >>> fitter = MultiNest(likelihood, priors)
    >>> result = fitter.run(seed=1)
    >>> result
    <MultiNest Result (median)
    teff_0    5718.882155
    logg_0       4.294858
    mh_0        -0.312654
    vrad_1      25.452319
    dtype: float64>

Keyword arguments to ``run`` are passed to ``pymultinest.run``. ``seed``
sets the seed of the MultiNest random number generator; without it the seed
is taken from the clock, and the result changes slightly from run to run. The
output above is that of the first fit in a new Python session; it made 9527
likelihood evaluations. ``n_live_points`` sets the number of live points,
which is 400 by default.

By default, MultiNest writes its files to a temporary directory, which is
deleted after the run. To keep them, give a directory with
``MultiNest(likelihood, priors, run_dir='fit1')``; the files are named
``specgrid_multinest_*`` unless ``prefix`` is given. The MultiNest library
from conda-forge limits the length of these paths; see :doc:`installation`.

The result
==========

``fitter.run`` returns a ``MultiNestResult``:

``result.median``, ``result.mean``
    Median and mean of the weighted posterior samples, as a pandas Series.
``result.calculate_sigmas(n)``
    For each parameter, the values that enclose the central part of the
    posterior that a Gaussian has within :math:`n\sigma`; for ``n=1``, the
    15.9 and 84.1 percentiles.
``result.maximum``
    The posterior sample with the highest likelihood.
``result.posterior_data``
    The weighted posterior samples as a pandas DataFrame, with the columns
    ``weights``, ``loglikelihood`` and one column per free parameter. The
    ``loglikelihood`` column holds :math:`-2 \ln L`, as written by MultiNest,
    not :math:`\ln L`.
``result.plot_triangle()``
    Corner plot of the posterior. It requires the corner package, which is
    in the ``multinest`` extra.

For the synthetic spectrum above, the 1σ intervals from
``result.calculate_sigmas(1)`` are 5644 to 5788 K for ``teff_0``, 4.16 to 4.43
for ``logg_0``, -0.371 to -0.260 for ``mh_0`` and 23.2 to 27.5 km/s for
``vrad_1``; each contains the input value.

``run`` does not change the parameters of the likelihood. To evaluate the
model at the posterior medians, set them::

    >>> for name, value in result.median.items():
    ...     setattr(likelihood, name, value)
    >>> model_wavelength, model_flux = model()

The model and the likelihood share their components, so setting the
parameters on ``likelihood`` also sets them on ``model``.

To save a result and read it back::

    >>> result.to_hdf('fit_result.h5')
    >>> from starkit.fitkit.multinest import MultiNestResult
    >>> result = MultiNestResult.from_hdf5('fit_result.h5')

``MultiNestResult.from_multinest_basename(basename, parameter_names)`` reads
the result from the files of a run kept in ``run_dir``, with ``basename`` the
directory and prefix, such as ``'fit1/specgrid_multinest'``.

Photometry and a spectrum together
==================================

``assemble_model`` with both ``spectrum`` and ``filter_set`` returns a model
whose outputs are the wavelength and flux of the model spectrum and the
magnitudes. A likelihood operation for each part is joined with ``&``, and
``Addition`` adds the two logarithms of the likelihood::

    from starkit.fitkit.likelihoods import Addition, PhotometryColorLikelihood

    sdss = ['SLOAN/SDSS.u', 'SLOAN/SDSS.g', 'SLOAN/SDSS.r', 'SLOAN/SDSS.i']
    model = assemble_model(grid, spectrum=spectrum, filter_set=sdss,
                           mag_type='ab', vrad=0., a_v=0., normalize_npol=3)
    likelihood = model | (SpectralChi2Likelihood(spectrum)
                          & PhotometryColorLikelihood(magnitudes,
                                                      uncertainties)) | Addition()

``magnitudes`` and ``uncertainties`` are the observed magnitudes and their
uncertainties, in the order of the filters. The extinction ``a_v`` is applied
to the magnitudes only. The example notebook
:doc:`examples/photometry_spectroscopy_fit` fits a synthetic star with
:math:`A_V` = 0.5 mag with its spectrum, its SDSS colors, and both: the
colors alone give :math:`A_V` = 0.81 +0.35 -0.36 mag, and the spectrum and
colors together give 0.49 +0.06 -0.06 mag.

MultiNest
=========

When using MultiNest please make sure to cite their works (arXiv:0809.3437,
arXiv:0704.3704, arXiv:1306.2144) as well as the wrapper PyMultiNest
(`Buchner et al. 2014, A&A <http://www.aanda.org/articles/aa/abs/2014/04/aa22971-13/aa22971-13.html>`_).
