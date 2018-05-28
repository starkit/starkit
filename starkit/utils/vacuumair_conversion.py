## Adding the vacuum to air conversion
# from http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion
from astropy import units as u


def convert_vacuum2air(wavelength):
    """
    Convert vacuum to air (using the formula given in
    http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion)
    Parameters
    ----------
    wavelength: float or astropy.Quantity
        if float is given wavelength is assumed to be in Angstrom

    Returns
    -------
        : astropy.Quantity
        wavelength in air
    """
    if hasattr(wavelength, 'unit'):
        wave_unit = wavelength.unit
    else:
        wave_unit = u.angstrom

    vac_wavelength = u.Quantity(wavelength, 'angstrom')
    s2 = 1e4 / vac_wavelength.value
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
    air_wavelength = (vac_wavelength / n).to(wave_unit)
    return air_wavelength

def convert_air2vacuum(wavelength):
    """
    Convert air to vacuum (using the formula given in
    http://www.astro.uu.se/valdwiki/Air-to-vacuum%20conversion)
    Parameters
    ----------
    wavelength: float or astropy.Quantity
        if float is given wavelength is assumed to be in Angstrom

    Returns
    -------
        : astropy.Quantity
        wavelength in vacuum
    """

    if hasattr(wavelength, 'unit'):
        wave_unit = wavelength.unit
    else:
        wave_unit = u.angstrom

    air_wavelength = u.Quantity(wavelength, 'angstrom')
    s2 = 1e4 / air_wavelength.value
    n = (1 + 0.00008336624212083 + 0.02408926869968 / (130.1065924522 - s2) +
         0.0001599740894897 / (38.92568793293 - s2))
    air_wavelength = (air_wavelength * n).to(wave_unit)
    return air_wavelength


