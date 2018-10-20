import numpy as np
from scipy import ndimage as nd
from scipy.interpolate import interp1d


from astropy import units as u
from astropy.io import fits

from starkit.gridkit.io.process import BaseProcessGrid
from starkit.gridkit.util import convolve_to_resolution


class PhoenixProcessGrid(BaseProcessGrid):

    uv_wavelength = (500, 3000)
    uv_R = 3000 / 0.1
    oir_wavelength = (3000, 25000)
    oir_R = 5e5
    nir_wavelength = (25000, 55000)
    nir_R = 1e5

    def __init__(self, index, input_wavelength, meta, wavelength_start=0*u.angstrom, wavelength_stop=np.inf*u.angstrom,
                 R=5000.0, R_sampling=4, pre_sampling=2):
        """
        Parameters
        ----------
        index: pandas.DataFrame
        input_wavelength : astropy.units.Quantity
        meta : pandas.Series
        wavelength_start : astropy.units.Quantity
        wavelength_stop : astropy.units.Quantity
        R : float
        R_sampling : integer
        pre_sampling: integer
            Select the sampling that you want to have before convolution (there are strange jumps in the diff)
        """
        super(PhoenixProcessGrid, self).__init__(index, input_wavelength, meta, wavelength_start=wavelength_start,
                                                 wavelength_stop=wavelength_stop, R=R, R_sampling=R_sampling)

        self.wavelength_interp = []
        for wl_region in ['uv', 'oir', 'nir']:
            current_wl = getattr(self, '{0}_wavelength'.format(wl_region))
            current_R = getattr(self, '{0}_R'.format(wl_region))
            wavelength_interp = self._logrange(current_wl[0], current_wl[1],
                                               current_R, pre_sampling)
            setattr(self, '{0}_wavelength_interp'.format(wl_region),
                    wavelength_interp)
            self.wavelength_interp.append(wavelength_interp)


        self.wavelength_interp = np.hstack(tuple(self.wavelength_interp))


    def interp_wavelength(self, flux):
        cut_flux = flux[self.start_idx:self.stop_idx]
        processed_wavelengths = []
        processed_fluxes = []
        for wl_region in ['uv', 'oir', 'nir']:
            current_R = getattr(self, '{0}_R'.format(wl_region))
            current_wl_interp = getattr(
                self, '{0}_wavelength_interp'.format(wl_region))
            rescaled_R = 1 / np.sqrt((1/self.R)**2 - (1/current_R)**2 )
            sigma = ((current_R / rescaled_R) * self.R_sampling /
                     (2 * np.sqrt(2 * np.log(2))))

            #The sampling is very strange (it's not really log - it's linear with jumps (every 5000 angstrom)
            #Specifically the UV is sampled at 0.1 angstrom which is also stated as the delta_lambda
            interp_flux = interp1d(self.cut_wavelength, cut_flux,
                                   bounds_error=False)(
                current_wl_interp
            )
            nan_filter = ~np.isnan(interp_flux)

            if np.all(~nan_filter):
                processed_wavelengths.append(np.array([]))
                processed_fluxes.append(np.array([]))
            else:
                processed_wavelengths.append(current_wl_interp[nan_filter])
                processed_fluxes.append(nd.gaussian_filter1d(
                    interp_flux[nan_filter], sigma))


        processed_wavelength = np.hstack(tuple(processed_wavelengths))
        processed_flux = np.hstack(tuple(processed_fluxes))

        if self.output_wavelength is None:
            initial_output_flux = interp1d(processed_wavelength, processed_flux,
                     bounds_error=False)(self.initial_output_wavelength)
            self.output_wavelength = self.initial_output_wavelength[~np.isnan(
                initial_output_flux)]
            output_flux = initial_output_flux[~np.isnan(initial_output_flux)]
        else:
            output_flux = interp1d(processed_wavelength, processed_flux)(
                self.output_wavelength
            )
        return output_flux

    @staticmethod
    def load_flux(fname):
        """
        load the 1D flux array from file

        Parameters
        ----------
        fname : str

        Returns
        -------
            : numpy.ndarray
        """

        #surface = fits.getval(fname, 'PHXREFF') ** 2 * 4 * np.pi
        flux = fits.getdata(fname).astype(np.float64)
        flux *= 1e-8 # converting from erg/s/cm to erg/s/angstrom
        return flux

