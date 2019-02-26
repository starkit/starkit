from starkit.utils import spectral

def test_fwhm2sigma():
    assert spectral.fwhm2sigma(spectral.FWHM2SIGMA_CONST) == 1.0
    assert spectral.sigma2fwhm(1.0) == spectral.FWHM2SIGMA_CONST 