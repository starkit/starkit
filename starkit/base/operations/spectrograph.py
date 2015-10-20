import warnings

import numpy as np
from numpy.polynomial import Polynomial
from scipy import ndimage as nd
from astropy import modeling
from astropy import units as u

from starkit.base.operations.base import (SpectralOperationModel,
                                          InstrumentOperationModel)
from starkit.fix_spectrum1d import Spectrum1D

from starkit.utils.spectral import prepare_observed

__all__ = ['InstrumentConvolveGrating', 'Interpolate', 'Normalize']

class SpectrographOperationModel(InstrumentOperationModel):
    pass

class InstrumentConvolveGrating(SpectrographOperationModel):
    """
    Convolve with a gaussian with given resolution to mimick an instrument
    assuming lambda / delta_lambda being constant

    Parameters
    ----------

    R : float or astropy.units.Quantity (unitless)
        resolution of the spectrum R = lambda/delta_lambda

    sampling: float
        number of pixels per resolution element (default=2.)

    """

    operation_name = 'resolution'

    R = modeling.Parameter()
    requires_observed_spectrum = False


    @classmethod
    def from_grid(cls, grid, R=np.inf):
        grid_R = getattr(grid, 'R', None)
        grid_sampling = getattr(grid, 'R_sampling', None)
        return cls(R=R, grid_R=grid_R, grid_sampling=grid_sampling)

    def __init__(self, R=np.inf, grid_R=None, grid_sampling=None):
        super(InstrumentConvolveGrating, self).__init__(R=R)
        self.grid_sampling = grid_sampling
        self.grid_R = grid_R

    def evaluate(self, wavelength, flux, R):
        if np.isinf(R):
            return wavelength, flux

        if self.grid_R is None:
            raise NotImplementedError('grid_R not given - this mode is not '
                                      'implemented yet')
        rescaled_R = 1 / np.sqrt((1/R)**2 - (1 / self.grid_R)**2 )

        sigma = ((self.grid_R / rescaled_R) * self.grid_sampling /
                     (2 * np.sqrt(2 * np.log(2))))

        return wavelength, nd.gaussian_filter1d(flux, sigma)

class InstrumentConvolveGrism(SpectrographOperationModel):
    """
    Convolve with a gaussian with given resolution to mimick an instrument
    assuming delta_lambda being constant

    Parameters
    ----------

    R : float or astropy.units.Quantity (unitless)
        resolution of the spectrum R = lambda/delta_lambda at wavelength

    sampling: float
        number of pixels per resolution element (default=2.)

    """

    operation_name = 'resolution'

    R = modeling.Parameter()
    requires_observed_spectrum = False

    @classmethod
    def from_grid(cls, wavelength, grid, R=np.inf):
        grid_R = getattr(grid, 'R', None)
        return cls(wavelength, R=R, grid_R=grid_R)


    def __init__(self, wavelength, R, grid_R, sampling=4, ):
        super(InstrumentConvolveGrism, self).__init__(R=R)
        self.wavelength = wavelength
        self.sampling = sampling
        self.fwhm2sigma = 1 / (2 * np.sqrt(np.log(2) * 2))
        self.grid_R = grid_R

    def evaluate(self, wavelength, flux, R):

        if np.isinf(R):
            return wavelength, flux

        rescaled_R = 1 / np.sqrt((1/R)**2 - (1 / self.grid_R)**2 )

        delta_lambda = (self.wavelength / rescaled_R) * self.fwhm2sigma


        new_wavelength = np.arange(wavelength[0], wavelength[-1],
                                   delta_lambda / self.sampling)

        new_flux = np.interp(new_wavelength, wavelength, flux)

        return new_wavelength, nd.gaussian_filter1d(new_flux, self.sampling)

class Interpolate(SpectrographOperationModel):

    """
    This class can be called to do a interpolation on a given spectrum.
    You must initialize it with the observed spectrum. The output will be a
    Spectrum1D object.

    Parameters
    ----------
    observed: Spectrum1D object
        This is the observed spectrum which you want to interpolate your
        (model) spectrum to.
    """

    requires_observed_spectrum = True
    operation_name = 'interpolate'

    def __init__(self, observed):
        super(SpectralOperationModel, self).__init__()
        self._update_observed_spectrum(observed)

    def _update_observed_spectrum(self, observed):
        self.observed = prepare_observed(observed)

    def evaluate(self, wavelength, flux):
        return self.observed.wavelength.value, np.interp(
            self.observed.wavelength.value, wavelength, flux)


class Normalize(SpectrographOperationModel):
    """Normalize a model spectrum to an observed one using a polynomial

    Parameters
    ----------
    observed : Spectrum1D object
        The observed spectrum to which the model should be matched
    npol : int
        The degree of the polynomial
    """

    requires_observed_spectrum = True

    operation_name = 'normalize'

    def __init__(self, observed, npol):
        super(Normalize, self).__init__()
        self.npol = npol
        self._update_observed_spectrum(observed)

    def _update_observed_spectrum(self, observed):
        self.observed = prepare_observed(observed)

        self.signal_to_noise = (self.observed.flux.value /
                                self.observed.uncertainty.value)
        self.flux_unit = observed.unit
        self._rcond = (len(observed.flux.value) *
                       np.finfo(observed.flux.dtype).eps)
        self._Vp = np.polynomial.polynomial.polyvander(
            observed.wavelength.value/observed.wavelength.mean().value - 1., self.npol)
        self.domain = u.Quantity([observed.wavelength.min().value,
                                  observed.wavelength.max().value])
        self.window = self.domain/observed.wavelength.mean().value - 1.


    def evaluate(self, wavelength, flux):
        # V[:,0]=mfi/e, Vp[:,1]=mfi/e*w, .., Vp[:,npol]=mfi/e*w**npol

        V = self._Vp * (flux / self.observed.uncertainty.value)[:, np.newaxis]
        # normalizes different powers
        scl = np.sqrt((V*V).sum(0))
        if np.isfinite(scl[0]):  # check for validity before evaluating
            sol, resids, rank, s = np.linalg.lstsq(V/scl, self.signal_to_noise,
                                                   self._rcond)
            sol = (sol.T / scl).T
            if rank != self._Vp.shape[-1] - 1:
                msg = "The fit may be poorly conditioned"
                warnings.warn(msg)

            fit = np.dot(V, sol) * self.observed.uncertainty.value
            # keep coefficients in case the outside wants to look at it
            self.polynomial = Polynomial(sol, domain=self.domain.value,
                                         window=self.window.value)
            return wavelength, fit
        else:
            return wavelength, flux


class NormalizeParts(SpectrographOperationModel):
    """Normalize a model spectrum to an observed one in multiple parts

    Here, different parts could, e.g., be different echelle orders or
    different chips for GMOS spectra.

    Parameters
    ----------
    observed : Spectrum1D object
        The observed spectrum to which the model should be matched
    parts : list of slices, index arrays, or boolean arrays
        Different parts of the observed spectrum which should be normalized
        separately.
    npol : list of int
        Polynomial degrees for the different parts
    """


    requires_observed_spectrum = True

    def __init__(self, observed, parts, npol):
        super(NormalizeParts, self).__init__()
        self.npol = npol

        self._update_observed_spectrum(observed, parts)

    def _update_observed_spectrum(self, observed_spectrum, parts):
        self.parts = parts
        self.normalizers = []
        try:
            if len(parts) != len(self.npol):
                raise ValueError("List of parts should match in length to "
                                 "list of degrees")
        except TypeError:  # npol is single number
            npol = [self.npol] * len(parts)
        else:
            npol = self.npol

        for part, _npol in zip(parts, npol):
            self.normalizers.append(
                Normalize(self.spectrum_1d_getitem(observed_spectrum, part),
                          _npol))


    @staticmethod
    def spectrum_1d_getitem(observed, part):
        observed_part = Spectrum1D.from_array(
            observed.wavelength[part],
            observed.flux[part])
        if getattr(observed, 'uncertainty', None) is not None:
            observed_part.uncertainty = getattr(observed.uncertainty, 'array',
                                                observed.uncertainty)[part]
        return observed_part


    def evaluate(self, wavelength, flux):
        fit_flux = np.zeros_like(flux)
        for part, normalizer in zip(self.parts, self.normalizers):
            _wave, fit_flux[part] = normalizer(wavelength[part], flux[part])

        return wavelength, fit_flux
