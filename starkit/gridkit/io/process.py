import numpy as np
from scipy.interpolate import interp1d
from scipy import ndimage as nd
from astropy import units as u

class BaseProcessGrid(object):

    def __init__(self, index, input_wavelength, meta, wavelength_start=0*u.angstrom, wavelength_stop=np.inf*u.angstrom,
                 R=None, R_sampling=4):

        self.index = index
        self.meta = meta
        self.R = R
        self.R_sampling=R_sampling

        self.wavelength_start = wavelength_start
        self.wavelength_stop = wavelength_stop

        input_wavelength = input_wavelength.to('angstrom')

        self.start_idx, self.stop_idx = self._get_wavelength_bounds(
            input_wavelength, wavelength_start, wavelength_stop)
        self.cut_wavelength = input_wavelength[
                              self.start_idx:self.stop_idx]

        self.output_wavelength = self._logrange(
            wavelength_start, wavelength_stop, R, R_sampling)






    @staticmethod
    def _get_wavelength_bounds(wavelength, wavelength_start, wavelength_stop):
        start_idx = np.max([wavelength.searchsorted(wavelength_start) - 1, 0])
        stop_idx = np.min([wavelength.searchsorted(wavelength_stop) + 1, len(wavelength)])
        return start_idx, stop_idx

    @staticmethod
    def _logrange(wavelength_start, wavelength_stop, R, sampling):
        return np.exp(
            np.arange(np.log(u.Quantity(wavelength_start, 'angstrom').value),
                      np.log(u.Quantity(wavelength_stop, 'angstrom').value),
                      1 / (sampling * float(R))))

