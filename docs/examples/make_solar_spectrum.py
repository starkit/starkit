"""Make starkit/data/sun_kurucz2005_r1500.fits, the solar spectrum used by
the example solar_spectrum_fit.ipynb.

The source is the solar irradiance atlas of Kurucz (2005), described in
Kurucz, R. L. 2006, arXiv:astro-ph/0605029: 300 to 1000 nm sampled every
0.025 nm, broadened with a Gaussian of 0.05 nm FWHM, in W m-2 nm-1, calibrated
to Thuillier et al. (2004), with the telluric lines divided out and the
gravitational redshift removed.

The wavelengths of the atlas file are vacuum wavelengths. The atlas poster
says the scans were joined in air, but the file was made from the wavenumber
file, and the cores of Na D1, Na D2, H alpha, H beta and Mg b 5183 lie within
1.6 km/s of their vacuum wavelengths and 84 km/s from their air wavelengths.

The atlas is smoothed to R = 1500, the resolving power of the test grid
phoenix_t4000_10000_w3000_9000_r3000.h5 (labeled R = 3000; see BUGS.md), and
sampled on the wavelengths of that grid. R is lambda / FWHM of a Gaussian
line-spread function.

Usage:
    python make_solar_spectrum.py [--input irradthuwl.br_05nmresamp] [--output file.fits]

Without --input, the atlas file is downloaded.
"""

import argparse
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.table import Table
from astropy.utils.data import download_file

SOURCE_URL = ('http://kurucz.harvard.edu/sun/irradiance2005/'
              'irradthuwl.br_05nmresamp')
REFERENCE = 'Kurucz, R. L. 2006, arXiv:astro-ph/0605029'
N_ATLAS_POINTS = 28001
ATLAS_FWHM = 0.5  # Angstrom, the Gaussian broadening of the atlas file
RESOLUTION = 1500.
# The wavelengths of the test grid: 3000 to 9000 Angstrom, R = 3000 with four
# points per resolution element (BaseProcessGrid._logrange).
GRID_START, GRID_STOP, GRID_R, GRID_SAMPLING = 3000., 9000., 3000., 4
FWHM_TO_SIGMA = 1 / (2 * np.sqrt(2 * np.log(2)))
DEFAULT_OUTPUT = (Path(__file__).resolve().parents[2] / 'starkit' / 'data'
                  / 'sun_kurucz2005_r1500.fits')


def read_atlas(path):
    """Wavelength (Angstrom) and irradiance (W m-2 nm-1) from the atlas file.

    The header lines and the line giving the number of points and the
    maximum have more or fewer than two numbers and are skipped.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 2:
                continue
            try:
                rows.append((float(parts[0]), float(parts[1])))
            except ValueError:
                continue
    wavelength_nm, flux = np.array(rows).T
    if wavelength_nm.size != N_ATLAS_POINTS:
        raise ValueError(f'expected {N_ATLAS_POINTS} points in {path}, '
                         f'found {wavelength_nm.size}')
    return wavelength_nm * 10., flux


def smooth_to_resolution(wavelength, flux, output_wavelength, resolution,
                         input_fwhm):
    """Smooth a spectrum to a resolving power and sample it.

    At each output wavelength, the result is the mean of the input points
    within five standard deviations, weighted with a Gaussian of FWHM
    sqrt((lambda / resolution)**2 - input_fwhm**2), so that a line of FWHM
    input_fwhm has FWHM lambda / resolution afterward.

    Parameters
    ----------
    wavelength, flux : numpy.ndarray
        the input spectrum, wavelength increasing
    output_wavelength : numpy.ndarray
        wavelengths at which the smoothed spectrum is sampled, in the unit of
        wavelength
    resolution : float
        resolving power lambda / FWHM
    input_fwhm : float
        FWHM of the Gaussian broadening of the input, in the unit of
        wavelength

    Returns
    -------
    numpy.ndarray
        the smoothed flux at output_wavelength
    """
    sigma = (np.sqrt((output_wavelength / resolution)**2 - input_fwhm**2)
             * FWHM_TO_SIGMA)
    low = np.searchsorted(wavelength, output_wavelength - 5 * sigma)
    high = np.searchsorted(wavelength, output_wavelength + 5 * sigma)
    smoothed = np.empty_like(output_wavelength)
    for i, (center, width, start, stop) in enumerate(
            zip(output_wavelength, sigma, low, high)):
        weights = np.exp(-0.5 * ((wavelength[start:stop] - center) / width)**2)
        smoothed[i] = np.sum(weights * flux[start:stop]) / np.sum(weights)
    return smoothed


def check_smoothing():
    """Smooth lines of FWHM ATLAS_FWHM on the atlas sampling and check that
    their FWHM becomes lambda / RESOLUTION to 1 %."""
    wavelength = np.arange(3000., 10000.0001, 0.25)
    for center in [4000., 6000., 8500.]:
        depth = np.exp(-0.5 * ((wavelength - center)
                               / (ATLAS_FWHM * FWHM_TO_SIGMA))**2)
        expected = center / RESOLUTION
        output = np.arange(center - 5 * expected, center + 5 * expected, 0.01)
        smoothed = smooth_to_resolution(wavelength, 1 - 0.5 * depth, output,
                                        RESOLUTION, ATLAS_FWHM)
        line = 1 - smoothed
        mean = np.sum(output * line) / np.sum(line)
        sigma = np.sqrt(np.sum((output - mean)**2 * line) / np.sum(line))
        fwhm = sigma / FWHM_TO_SIGMA
        if abs(fwhm / expected - 1) > 0.01:
            raise RuntimeError(f'line at {center} A: FWHM {fwhm:.4f} A after '
                               f'smoothing, expected {expected:.4f} A')
        print(f'line at {center:.0f} A: FWHM {fwhm:.4f} A, expected '
              f'{expected:.4f} A')


def grid_wavelength():
    return np.exp(np.arange(np.log(GRID_START), np.log(GRID_STOP),
                            1 / (GRID_SAMPLING * GRID_R)))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--input', help='the atlas file; downloaded if not '
                        'given')
    parser.add_argument('--output', default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()

    check_smoothing()
    path = args.input or download_file(SOURCE_URL, cache=False)
    wavelength, flux = read_atlas(path)

    output_wavelength = grid_wavelength()
    # Keep the output points whose smoothing kernel lies inside the atlas.
    output_wavelength = output_wavelength[
        (output_wavelength > wavelength[0] + 10.)
        & (output_wavelength < wavelength[-1] - 10.)]
    output_flux = smooth_to_resolution(wavelength, flux, output_wavelength,
                                       RESOLUTION, ATLAS_FWHM)

    table = Table([output_wavelength * u.angstrom,
                   output_flux * u.W / u.m**2 / u.nm],
                  names=['WAVELENGTH', 'FLUX'])
    table.meta['OBJECT'] = 'Sun'
    table.meta['REFERENC'] = REFERENCE
    table.meta['SOURCE'] = SOURCE_URL
    table.meta['WAVETYPE'] = 'vacuum'
    table.meta['RESOLUTN'] = RESOLUTION
    table.meta['HISTORY'] = [
        'Solar irradiance atlas of Kurucz (2005), 300-1000 nm, calibrated to',
        'Thuillier et al. (2004), telluric lines divided out, gravitational',
        'redshift removed. Vacuum wavelengths (line cores within 1.6 km/s of',
        'their vacuum wavelengths). Smoothed with a Gaussian to R = 1500',
        '(R = lambda / FWHM), taking the 0.05 nm FWHM broadening of the atlas',
        'file into account, and sampled on the wavelengths of the test grid',
        'phoenix_t4000_10000_w3000_9000_r3000.h5 (labeled R = 3000, effective',
        'R = 1500). Made with docs/examples/make_solar_spectrum.py.',
    ]
    table.write(args.output, overwrite=True)
    print(f'{len(table)} points from {output_wavelength[0]:.2f} to '
          f'{output_wavelength[-1]:.2f} Angstrom written to {args.output}')


if __name__ == '__main__':
    main()
