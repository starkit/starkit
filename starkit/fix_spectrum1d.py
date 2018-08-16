import numpy as np
from astropy import units as u

class SKSpectrum1D(object):

    @classmethod
    def from_array(cls, wavelength, flux, uncertainty=None):
        """
        Parameters
        ----------
        wavelength: astropy.units.Quantity
        flux: astropy quantity
        uncertainty: astropy quantity

        Returns
        -------
            : SKSpectrum1D object
        """

        return cls(wavelength, flux, uncertainty)

    def __init__(self, wavelength, flux, uncertainty):
        self.wavelength = wavelength
        self.flux = flux
        self.uncertainty = uncertainty


    def uncertainty_getter(self):
        return self._uncertainty

    def uncertainty_setter(self, value):
        if value is None:
            self._uncertainty = None
        else:
            self._uncertainty = u.Quantity(value, self.flux.unit)

    uncertainty = property(uncertainty_getter, uncertainty_setter)

    def slice_index(self, start=None, stop=None, step=None):
        """
        Slice the spectrum according to index

        Parameters
        ----------
        start: int
        stop: int
        step: int

        Returns
        -------
            : SKSpectrum
        """

        spectral_slice = slice(start, stop, step)
        if self.uncertainty is not None:
            new_uncertainty = self.uncertainty[spectral_slice]
        else:
            new_uncertainty = None

        return SKSpectrum1D(self.wavelength[spectral_slice], self.flux[spectral_slice],
                            new_uncertainty)

    def slice_wavelength(self, start_wavelength=None, stop_wavelength=None):
        """
        Slice the spectrum according to wavelength

        Parameters
        ----------
        start_wavelength: astropy.units.Quantity or float
            if float is given the assumed wavelength unit will be that of the
            spectrum
        stop_wavelength
            if float is given the assumed wavelength unit will be that of the
            spectrum

        Returns
        -------
            : SKSpectrum
        """
        if start_wavelength is None:
            start = None
        else:
            start_wavelength = u.Quantity(start_wavelength, unit=self.wavelength.unit)
            start = self.wavelength.searchsorted(start_wavelength, side='left')

        if stop_wavelength is None:
            stop = None
        else:
            stop_wavelength = u.Quantity(stop_wavelength, unit=self.wavelength.unit)
            stop = self.wavelength.searchsorted(stop_wavelength, side='right')

        return self.slice_index(start=start, stop=stop)

    def get_nan_cleaned(self):
        """
        Get a new spectrum that has all NaNs removed

        Returns
        -------
            : SKSpectrum1D

        """

        nan_filter = ~np.isnan(self.wavelength)
        nan_filter = nan_filter & ~np.isnan(self.flux)
        if self.uncertainty is not None:
            nan_filter = nan_filter & ~np.isnan(self.uncertainty)
            nan_filter = nan_filter & (self.uncertainty > 0.0)

        if self.uncertainty is None:
            new_uncertainty = None
        else:
            new_uncertainty = self.uncertainty[nan_filter]

        return SKSpectrum1D(self.wavelength[nan_filter], self.flux[nan_filter],
                            new_uncertainty)