from astropy import constants as const, units as u
from astropy import modeling
import scipy.ndimage as nd

__all__ = ['RotationalBroadening', 'CCM89Extinction', 'DopplerShift']

import numpy as np

from starkit.base.operations.base import SpectralOperationModel

class StellarOperationModel(SpectralOperationModel):
    pass

class RotationalBroadening(StellarOperationModel):
	"""
	The rotational broadening kernel was taken from 
	Observation and Analysis of Stellar Photospheres
	by David Gray
	"""
    operation_name = 'rotation'
    vrot = modeling.Parameter()
    limb_darkening = modeling.Parameter(fixed=True, default=0.6)

    @classmethod
    def from_grid(cls, grid, vrot=0):
        velocity_per_pix = getattr(grid, 'velocity_per_pix', None)

        return cls(velocity_per_pix=velocity_per_pix, vrot=vrot)



    def __init__(self, velocity_per_pix=None, vrot=0):
        super(RotationalBroadening, self).__init__(vrot=vrot)

        self.c_in_kms = const.c.to(u.km / u.s).value

        if velocity_per_pix is not None:
            self.log_sampling = True
            self.velocity_per_pix = u.Quantity(velocity_per_pix, u.km / u.s).value
        else:
            self.log_sampling = False
            self.velocity_per_pix = None

    def rotational_profile(self, vrot, limb_darkening):
        vrot = float(vrot)
        limb_darkening = float(limb_darkening)
        vrot_by_c = np.maximum(0.0001, np.abs(vrot)) / self.c_in_kms

        half_width_pix = np.round((vrot /
                                   self.velocity_per_pix)).astype(int)
        profile_velocity = (np.linspace(-half_width_pix, half_width_pix,
                                       2 * half_width_pix + 1)
                            * self.velocity_per_pix)
        profile = np.maximum(0.,
                             1. - (profile_velocity / vrot) ** 2)

        profile = ((2 * (1 - limb_darkening) * np.sqrt(profile) +
                    0.5 * np.pi * limb_darkening * profile) /
                   (np.pi * vrot_by_c * (1. - limb_darkening / 3.)))
        return profile / profile.sum()


    def evaluate(self, wavelength, flux, v_rot, limb_darkening):
        v_rot = np.asscalar(v_rot)
        limb_darkening = np.asscalar(limb_darkening)

        if self.velocity_per_pix is None:
            raise NotImplementedError('Regridding not implemented yet')

        if np.abs(v_rot) < 1e-5:
            return wavelength, flux

        profile = self.rotational_profile(v_rot, limb_darkening)

        return wavelength, nd.convolve1d(flux, profile)

class DopplerShift(StellarOperationModel):

    operation_name = 'doppler'

    vrad = modeling.Parameter()

    def __init__(self, vrad):
        super(DopplerShift, self).__init__(vrad=vrad)
        self.c_in_kms = const.c.to(u.km / u.s).value


    def evaluate(self, wavelength, flux, vrad):
        beta = vrad / self.c_in_kms
        doppler_factor = np.sqrt((1+beta) / (1-beta))
        return wavelength * doppler_factor, flux



class CCM89Extinction(StellarOperationModel):

    operation_name = 'ccm89_extinction'

    a_v = modeling.Parameter(default=0.0)
    r_v = modeling.Parameter(default=3.1, fixed=True)

    @property
    def ebv(self):
        return self.a_v / self.r_v

    def __init__(self, a_v=0.0, r_v=3.1):
        super(CCM89Extinction, self).__init__(a_v=a_v, r_v=r_v)


    def evaluate(self, wavelength, flux, a_v, r_v):
        from specutils import extinction
        wavelength = np.array(wavelength)
        extinction_factor = np.ones_like(wavelength)
        valid_wavelength = ((wavelength > 910) & (wavelength < 33333))

        extinction_factor[valid_wavelength] = 10 ** (-0.4 * extinction.extinction_ccm89(
            wavelength[valid_wavelength] * u.angstrom, a_v=np.abs(a_v),
            r_v=np.abs(r_v)))

        return wavelength, extinction_factor * flux

class Distance(StellarOperationModel):

    operation_name = 'distance'

    distance = modeling.Parameter(default=10.0) # in pc

    @classmethod
    def from_grid(cls, grid, distance=10.):
        lum_density2cgs = grid.flux_unit.to('erg / (s * angstrom)')
        return cls(distance=distance, lum_density2cgs=lum_density2cgs)

    def __init__(self, distance, lum_density2cgs=1.):
        super(Distance, self).__init__(distance=distance)
        self.pc2cm = u.pc.to(u.cm)
        self.lum_density2cgs = lum_density2cgs


    def evaluate(self, wavelength, flux, distance):
        conversion = self.lum_density2cgs / (4 * np.pi *
                                             (distance * self.pc2cm)**2)
        return wavelength, flux * conversion