Data directory
==============

This directory contains data files included with the affiliated package source
code distribution. Note that this is intended only for relatively small files
- large files should be externally hosted and downloaded as needed.

``alpha_lyr_stis_011.fits``
    The CALSPEC spectrum of Vega (Bohlin, Gordon & Tremblay 2014, PASP 126,
    711), from https://ssb.stsci.edu/trds/calspec/alpha_lyr_stis_011.fits,
    used for Vega zero points by ``starkit.filters`` (``VEGA_FILE``).

``sun_kurucz2005_r1500.fits``
    The solar irradiance atlas of Kurucz (2005; Kurucz, R. L. 2006,
    arXiv:astro-ph/0605029), from 3010 to 9000 Å in W m\ :sup:`-2` nm\ :sup:`-1`
    with vacuum wavelengths, smoothed to R = 1500 and sampled on the
    wavelengths of the documentation test grid. Used by the example
    ``docs/examples/solar_spectrum_fit.ipynb`` and made by
    ``docs/examples/make_solar_spectrum.py``, which records the source and the
    processing.
