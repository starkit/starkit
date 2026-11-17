"""Regression tests for the stellar operations.

astropy.modeling passes parameters to ``evaluate`` as length-1 arrays. Since
numpy 2, ``float()`` of such an array raises ``TypeError``, so every
operation is exercised through a model call rather than ``evaluate`` alone.
"""

from types import SimpleNamespace

import numpy as np
import pytest
from astropy import units as u

from starkit.base.operations.stellar import (
    CCM89Extinction,
    Distance,
    DopplerShift,
    RadialVelocity,
    RotationalBroadening,
)

C_KMS = 299792.458

# A_lambda for A_V = 1 from extinction.ccm89 (K. Barbary), the same CCM89
# algorithm that specutils.extinction.extinction_ccm89 used before starkit
# computed it internally.
CCM89_REFERENCE = [
    (22000.0, 3.1, 0.1135222130852181),   # IR, x = 0.455 um^-1
    (12000.0, 3.1, 0.3012309274775814),   # IR, x = 0.833 um^-1
    (5494.5, 3.1, 1.000001151914559),     # optical, x = 1.820 um^-1
    (4400.0, 3.1, 1.3245212536961641),    # optical, x = 2.273 um^-1
    (2175.0, 3.1, 3.1851012753144046),    # UV, x = 4.598 um^-1
    (1500.0, 3.1, 2.6638792147137162),    # UV with F_a, F_b, x = 6.667 um^-1
    (1100.0, 3.1, 4.2172065388623645),    # far UV, x = 9.091 um^-1
    (950.0, 3.1, 6.045613491104225),      # far-UV curve extrapolated past x = 10
    (22000.0, 5.0, 0.13167452735577523),
    (12000.0, 5.0, 0.3493980510296897),
    (5494.5, 5.0, 1.0000008365927981),
    (4400.0, 5.0, 1.2011076421249995),
    (2175.0, 5.0, 1.9742861111897607),
    (1500.0, 5.0, 1.4960346187487024),
    (1100.0, 5.0, 1.9740151765589788),
    (950.0, 5.0, 2.6410327744569178),
]


def _extinction_mag(wavelength, a_v, r_v):
    wavelength = np.atleast_1d(np.asarray(wavelength, dtype=float))
    model = CCM89Extinction(a_v=a_v, r_v=r_v)
    out_wavelength, out_flux = model(wavelength, np.ones_like(wavelength))
    np.testing.assert_array_equal(out_wavelength, wavelength)
    return -2.5 * np.log10(out_flux)


@pytest.mark.parametrize("wavelength, r_v, a_lambda", CCM89_REFERENCE)
def test_ccm89_matches_reference_curve(wavelength, r_v, a_lambda):
    np.testing.assert_allclose(_extinction_mag(wavelength, 1.0, r_v), [a_lambda], rtol=1e-10)


def test_ccm89_scales_linearly_with_a_v():
    np.testing.assert_allclose(_extinction_mag(22000.0, 30.0, 3.1), [30 * 0.1135222130852181], rtol=1e-10)


def test_ccm89_uses_absolute_values_of_a_v_and_r_v():
    wavelength = [4400.0, 22000.0]
    np.testing.assert_allclose(_extinction_mag(wavelength, -2.0, -3.1), _extinction_mag(wavelength, 2.0, 3.1))


@pytest.mark.parametrize("wavelength", [500.0, 900.0, 33333.2, 40000.0, 1e5])
def test_ccm89_leaves_flux_unchanged_outside_its_range(wavelength):
    # Extinction is applied between 910 and 33333 Angstrom, the range starkit
    # used with specutils.
    np.testing.assert_array_equal(_extinction_mag(wavelength, 5.0, 3.1), [0.0])


def test_doppler_shift_applies_relativistic_factor():
    wavelength = np.array([10000.0, 20000.0])
    flux = np.array([1.0, 2.0])
    out_wavelength, out_flux = DopplerShift(vrad=0.001 * C_KMS)(wavelength, flux)
    # sqrt((1 + beta) / (1 - beta)) for beta = 0.001
    np.testing.assert_allclose(out_wavelength, wavelength * 1.0010005005003753, rtol=1e-14)
    np.testing.assert_array_equal(out_flux, flux)


def test_radial_velocity_applies_classical_factor():
    wavelength = np.array([10000.0, 20000.0])
    out_wavelength, _ = RadialVelocity(vrad=0.001 * C_KMS)(wavelength, np.ones(2))
    np.testing.assert_allclose(out_wavelength, wavelength * 1.001, rtol=1e-14)


def test_distance_scales_flux_by_inverse_square_law():
    wavelength = np.array([10000.0, 20000.0])
    flux = np.array([1.0, 3.0])
    _, flux_10pc = Distance(distance=10.0)(wavelength, flux)
    _, flux_20pc = Distance(distance=20.0)(wavelength, flux)
    # 1 / (4 pi (10 pc in cm)^2)
    np.testing.assert_allclose(flux_10pc, flux * 8.357743635931359e-41, rtol=1e-12)
    np.testing.assert_allclose(flux_20pc, flux_10pc / 4.0, rtol=1e-14)


def test_distance_from_grid_converts_the_grid_flux_unit():
    # assemble_model builds Distance with from_grid, which reads the unit of
    # the grid fluxes. 1 W / nm is 1e6 erg / (s Angstrom).
    wavelength = np.array([10000.0, 20000.0])
    flux = np.array([1.0, 3.0])
    grid = SimpleNamespace(fluxes=np.ones((2, 2)) * u.Unit("W / nm"))
    _, flux_10pc = Distance.from_grid(grid, distance=10.0)(wavelength, flux)
    np.testing.assert_allclose(flux_10pc, flux * 1e6 * 8.357743635931359e-41, rtol=1e-12)


def test_operations_reject_multi_valued_parameters():
    with pytest.raises(ValueError, match="vrad"):
        DopplerShift(vrad=1.0).evaluate(np.ones(3), np.ones(3), np.array([1.0, 2.0]))


def _rotate_single_pixel_line(vrot):
    # 2 km/s per pixel; one absorption pixel of depth 0.5 in the middle.
    flux = np.ones(101)
    flux[50] = 0.5
    model = RotationalBroadening(velocity_per_pix=2.0, vrot=vrot)
    return model(np.arange(101.0), flux)


def test_rotational_broadening_spreads_line_over_kernel_width():
    _, flux = _rotate_single_pixel_line(20.0)
    # The kernel half-width is 20 / 2 = 10 pixels; its outermost samples are
    # zero, so 19 pixels are affected.
    absorbed = np.flatnonzero(flux < 1.0 - 1e-12)
    np.testing.assert_array_equal(absorbed, np.arange(41, 60))
    # Gray's profile with limb darkening 0.6, sampled at 21 points and
    # normalized: the central weight is 0.0698303146447221.
    np.testing.assert_allclose(1.0 - flux[50], 0.5 * 0.0698303146447221, rtol=1e-12)


def test_rotational_broadening_conserves_equivalent_width():
    _, flux = _rotate_single_pixel_line(20.0)
    np.testing.assert_allclose(np.sum(1.0 - flux), 0.5, rtol=1e-12)


def test_rotational_broadening_is_skipped_for_zero_velocity():
    _, flux = _rotate_single_pixel_line(0.0)
    assert flux[50] == 0.5
    np.testing.assert_array_equal(np.flatnonzero(flux < 1.0), [50])
