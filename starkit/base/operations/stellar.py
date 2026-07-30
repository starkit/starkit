from astropy import constants as const, units as u
from astropy import modeling
import scipy.ndimage as nd

__all__ = ['RotationalBroadening', 'CCM89Extinction', 'DopplerShift']

import numpy as np

def _as_scalar_float(x, name):
    arr = np.asarray(x)
    if arr.ndim == 0:
        return float(arr)
    if arr.size == 1:
        return float(arr.reshape(()))
    raise ValueError(
        f"{name} must be a scalar or length-1 array, got shape {arr.shape} and value {x!r}"
    )


def _ccm89_a_b(x):
    """
    Return the CCM89 a(x), b(x) coefficients for x in inverse microns.
    """
    x = np.asarray(x, dtype=float)
    a = np.zeros_like(x)
    b = np.zeros_like(x)

    ir = (x >= 0.3) & (x < 1.1)
    if np.any(ir):
        x_ir = x[ir]
        a[ir] = 0.574 * x_ir ** 1.61
        b[ir] = -0.527 * x_ir ** 1.61

    optical = (x >= 1.1) & (x < 3.3)
    if np.any(optical):
        y = x[optical] - 1.82
        a[optical] = (
            1.0
            + 0.17699 * y
            - 0.50447 * y ** 2
            - 0.02427 * y ** 3
            + 0.72085 * y ** 4
            + 0.01979 * y ** 5
            - 0.77530 * y ** 6
            + 0.32999 * y ** 7
        )
        b[optical] = (
            1.41338 * y
            + 2.28305 * y ** 2
            + 1.07233 * y ** 3
            - 5.38434 * y ** 4
            - 0.62251 * y ** 5
            + 5.30260 * y ** 6
            - 2.09002 * y ** 7
        )

    uv = (x >= 3.3) & (x <= 8.0)
    if np.any(uv):
        x_uv = x[uv]
        a_uv = 1.752 - 0.316 * x_uv - 0.104 / ((x_uv - 4.67) ** 2 + 0.341)
        b_uv = -3.090 + 1.825 * x_uv + 1.206 / ((x_uv - 4.62) ** 2 + 0.263)
        f_uv = x_uv >= 5.9
        if np.any(f_uv):
            y = x_uv[f_uv] - 5.9
            a_uv[f_uv] += -0.04473 * y ** 2 - 0.009779 * y ** 3
            b_uv[f_uv] += 0.2130 * y ** 2 + 0.1207 * y ** 3
        a[uv] = a_uv
        b[uv] = b_uv

    fuv = (x > 8.0) & (x <= 10.0)
    if np.any(fuv):
        y = x[fuv] - 8.0
        a[fuv] = -1.073 - 0.628 * y + 0.137 * y ** 2 - 0.070 * y ** 3
        b[fuv] = 13.670 + 4.257 * y - 0.420 * y ** 2 + 0.374 * y ** 3

    return a, b


def _ccm89_extinction_factor(wavelength, a_v, r_v):
    """
    Compute the CCM89 attenuation factor for wavelengths in Angstrom.
    """
    wavelength = np.asarray(wavelength, dtype=float)
    factor = np.ones_like(wavelength, dtype=float)
    valid = wavelength > 0.0
    if not np.any(valid):
        return factor

    x = 1.0 / (wavelength[valid] * 1e-4)  # inverse microns
    ccm_valid = (x >= 0.3) & (x <= 10.0)
    if not np.any(ccm_valid):
        return factor

    a, b = _ccm89_a_b(x[ccm_valid])
    a_lambda = np.abs(a_v) * (a + b / np.abs(r_v))
    valid_indices = np.flatnonzero(valid)
    factor[valid_indices[ccm_valid]] = 10.0 ** (-0.4 * a_lambda)
    return factor


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
        vrot = _as_scalar_float(vrot, "vrot")
        limb_darkening = _as_scalar_float(limb_darkening, "limb_darkening")
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
        v_rot = _as_scalar_float(v_rot, "v_rot")
        limb_darkening = _as_scalar_float(limb_darkening, "limb_darkening")

        if self.velocity_per_pix is None:
            raise NotImplementedError('Regridding not implemented yet')

        if np.abs(v_rot) < 1e-5:
            return wavelength, flux

        profile = self.rotational_profile(v_rot, limb_darkening)

        return wavelength, nd.convolve1d(flux, profile)

class DopplerShift(StellarOperationModel):
    '''
    Module to model the velocity of the star with the relativistic redshift correction. 
    '''

    operation_name = 'doppler'

    vrad = modeling.Parameter()

    def __init__(self, vrad):
        super(DopplerShift, self).__init__(vrad=vrad)
        self.c_in_kms = const.c.to(u.km / u.s).value


    def evaluate(self, wavelength, flux, vrad):
        vrad = _as_scalar_float(vrad, "vrad")
        beta = vrad / self.c_in_kms
        doppler_factor = np.sqrt((1+beta) / (1-beta))
        return wavelength * doppler_factor, flux

class RadialVelocity(StellarOperationModel):
    '''
    Module to model the classical definition of radial velocity
    '''
    operation_name = 'rv'

    vrad = modeling.Parameter()

    def __init__(self, vrad):
        super(RadialVelocity, self).__init__(vrad=vrad)
        self.c_in_kms = const.c.to(u.km / u.s).value


    def evaluate(self, wavelength, flux, vrad):
        vrad = _as_scalar_float(vrad, "vrad")
        beta = vrad / self.c_in_kms
        doppler_factor = (1.0 + beta)
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
        a_v = _as_scalar_float(a_v, "a_v")
        r_v = _as_scalar_float(r_v, "r_v")
        wavelength = np.asarray(wavelength, dtype=float)
        extinction_factor = _ccm89_extinction_factor(wavelength, a_v, r_v)
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
        distance = _as_scalar_float(distance, "distance")
        conversion = self.lum_density2cgs / (4 * np.pi *
                                             (distance * self.pc2cm)**2)
        return wavelength, flux * conversion
