from starkit.utils import spectral
from astropy.units import Quantity
from starkit.fix_spectrum1d import SKSpectrum1D


def test_fwhm2sigma():
    assert spectral.fwhm2sigma(spectral.FWHM2SIGMA_CONST) == 1.0
    assert spectral.sigma2fwhm(1.0) == spectral.FWHM2SIGMA_CONST


def test_prepare_observed():
    def check(spec, obj):
        if obj.uncertainty is not None:
            if (type(obj.wavelength.value) == float) | (type(obj.wavelength.value) == int):
                assert spec.wavelength.value == obj.wavelength.value
                assert spec.flux.value == obj.flux.value
                assert spec.uncertainty.value == obj.uncertainty.value
            else:
                for i in range(len(spec.wavelength.value)):
                    assert spec.wavelength.value[i] == obj.wavelength.value[i]
                    assert spec.flux.value[i] == obj.flux.value[i]
                    assert spec.uncertainty.value[i] == obj.uncertainty.value[i]
        else:
            if (type(obj.wavelength.value) == float) | (type(obj.wavelength.value) == int):
                assert spec.wavelength.value == obj.wavelength.value
                assert spec.flux.value == obj.flux.value
                assert spec.uncertainty is not None
            else:
                for i in range(len(spec.wavelength.value)):
                    assert spec.wavelength.value[i] == obj.wavelength.value[i]
                    assert spec.flux.value[i] == obj.flux.value[i]
                    assert spec.uncertainty is not None
        assert type(spec) == SKSpectrum1D

    w = Quantity(1000, 'Angstrom')
    f = Quantity(100, 'Jy')
    u = Quantity(3, 'Jy')
    obj = SKSpectrum1D(w, f, u)
    spec = spectral.prepare_observed(obj)
    check(spec, obj)

    w = Quantity([1000, 1500, 2000], 'Angstrom')
    f = Quantity([100, 150, 100], 'Jy')
    u = Quantity([2, 5, 1], 'Jy')
    obj = SKSpectrum1D(w, f, u)
    spec = spectral.prepare_observed(obj)
    check(spec, obj)

    w = Quantity([1000, 1500, 2000], 'Angstrom')
    f = Quantity([100, 150, 100], 'Jy')
    obj = SKSpectrum1D(w, f, None)
    spec = spectral.prepare_observed(obj)
    check(spec, obj)
