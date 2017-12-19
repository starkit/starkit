import numpy as np
from scipy.interpolate import interp1d
from scipy import ndimage as nd
from astropy import units as u

class BaseProcessGrid(object):

    def __init__(self, input_wavelength, wavelength_range = (0, np.inf), R=None, R_sampling=4):

        input_wavelength = input_wavelength.to('angstrom').value
        self.start_idx, self.stop_idx = self._get_wavelength_bounds(
            input_wavelength, wavelength_range)
        self.input_wavelength = input_wavelength
        output_wavelength = self._logrange(
            input_wavelength[self.start_idx],
            input_wavelength[self.stop_idx - 1], R, sampling)

        self.output_wavelength = self.test_interpolation(input_wavelength,
                                                         output_wavelength)

    def test_interpolation(self, input_wavelength, output_wavelength,
                           epsilon=1e-10):
        test_flux = np.ones_like(input_wavelength)

        interp_test_flux = self.process(test_flux, output_wavelength,
                                        bounds_error=False)

        if np.any(np.isnan(interp_test_flux)):
            interp_test_flux = self.process(test_flux, output_wavelength
                                            + epsilon, bounds_error=False)
            if not np.any(np.isnan(interp_test_flux)):
                output_wavelength = output_wavelength + epsilon

        return output_wavelength[~np.isnan(interp_test_flux)]





    @staticmethod
    def _get_wavelength_bounds(wavelength, wavelength_range):
        start_idx = wavelength.searchsorted(
            wavelength_range[0])
        stop_idx = wavelength.searchsorted(wavelength_range[1])
        return start_idx, stop_idx

    @staticmethod
    def _logrange(wavelength_start, wavelength_end, R, sampling):
        return np.exp(
            np.arange(np.log(wavelength_start), np.log(wavelength_end),
                      1 / (sampling * float(R))))


    def process(self, flux, output_wavelength=None, bounds_error=True):
        if output_wavelength is None:
            output_wavelength = self.output_wavelength
        return interp1d(self.input_wavelength, flux,
                        bounds_error=bounds_error)(output_wavelength)

