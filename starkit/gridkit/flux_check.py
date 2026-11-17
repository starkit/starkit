"""Compare the flux scale of a grid with the surface flux of a blackbody.

The ``flux_unit`` label of a grid does not always describe its fluxes: grids
written by StarKit before January 2025 are labeled erg/s/angstrom whether they
hold luminosity density or surface flux. ``check_flux_scale`` integrates a few
model spectra of a grid over its wavelength range and divides by the integral
of the surface flux of a blackbody at the same effective temperature,
pi * B_lambda(teff). The ratio identifies what the fluxes are:

- about 1 for surface flux, F_lambda (erg / (s cm2 Angstrom));
- about 1/pi for surface flux divided by pi;
- about 1/(4 pi) for the Eddington flux, F_lambda / (4 pi);
- 4 pi R**2, with R the stellar radius in cm (about 6e22 for the Sun), for
  luminosity density, L_lambda (erg / (s Angstrom));
- continuum-normalized spectra have values near 1 and carry no flux scale.

Model spectra differ from a blackbody by tens of percent, more so over a
narrow wavelength range, so the ratios of surface fluxes found in practice
lie between about 0.5 and 2.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as u
from scipy.integrate import trapezoid

from starkit.gridkit.base import _read_grid_meta

__all__ = ['ModelFluxCheck', 'FluxScaleCheck', 'blackbody_surface_flux',
           'check_flux_scale', 'format_flux_scale_check']

SURFACE_FLUX = u.erg / (u.s * u.cm ** 2 * u.AA)
LUMINOSITY_DENSITY = u.erg / (u.s * u.AA)

SURFACE = 'surface flux'
SURFACE_PI = 'surface flux divided by pi'
EDDINGTON = 'Eddington flux (surface flux divided by 4 pi)'
LUMINOSITY = 'luminosity density'
NORMALIZED = 'continuum-normalized'
UNRECOGNIZED = 'not recognized'
INVALID = 'invalid (non-finite or negative fluxes)'


@dataclass(frozen=True)
class ModelFluxCheck:
    """Flux scale of one model spectrum of a grid.

    ``ratio`` is the integral of the grid flux, in the unit of the label
    converted to erg/(s cm2 Angstrom) or erg/(s Angstrom), divided by the
    integral of pi * B_lambda(teff) over the same wavelengths. For
    luminosity density, ``implied_radius`` is the stellar radius in solar
    radii for which the ratio equals 4 pi R**2. A model with non-finite or
    negative fluxes is reported as invalid and its ratio is computed from
    its finite pixels only.
    """
    parameters: dict
    ratio: float
    median_flux: float
    interpretation: str
    implied_radius: float | None
    n_nonfinite: int = 0
    min_flux: float = 0.0


@dataclass(frozen=True)
class FluxScaleCheck:
    path: str
    flux_unit: str
    label_kind: str
    wavelength_min: float
    wavelength_max: float
    models: list
    interpretation: str
    consistent_with_label: bool | None


def blackbody_surface_flux(wavelength, teff):
    """Surface flux pi * B_lambda of a blackbody, in erg / (s cm2 Angstrom).

    Parameters
    ----------
    wavelength : numpy.ndarray
        Wavelengths in Angstrom.
    teff : float
        Temperature in K.
    """
    wavelength_cm = np.asarray(wavelength, dtype=float) * 1e-8
    h, c, k = const.h.cgs.value, const.c.cgs.value, const.k_B.cgs.value
    # The exponential overflows far on the Wien side, where B_lambda is 0.
    with np.errstate(over='ignore'):
        b_lambda = (2 * h * c ** 2 / wavelength_cm ** 5
                    / np.expm1(h * c / (wavelength_cm * k * teff)))
    return np.pi * b_lambda * 1e-8


def _label_kind(unit):
    if unit.physical_type == 'dimensionless':
        return 'dimensionless'
    if unit.is_equivalent(SURFACE_FLUX):
        return SURFACE
    if unit.is_equivalent(LUMINOSITY_DENSITY):
        return LUMINOSITY
    return 'other'


def _to_label_cgs(flux, unit, label_kind):
    if label_kind == SURFACE:
        return (flux * unit).to_value(SURFACE_FLUX)
    if label_kind == LUMINOSITY:
        return (flux * unit).to_value(LUMINOSITY_DENSITY)
    return flux * unit.scale if label_kind == 'dimensionless' else flux


def _interpret(ratio, median_flux, max_flux):
    if 0.3 <= median_flux <= 1.5 and max_flux <= 2.0:
        return NORMALIZED
    if ratio > 1e15:
        return LUMINOSITY
    if 0.5 <= ratio <= 2.0:
        return SURFACE
    if 0.2 <= ratio < 0.5:
        return SURFACE_PI
    if 0.05 <= ratio < 0.2:
        return EDDINGTON
    return UNRECOGNIZED


_CONSISTENT = {SURFACE: {SURFACE}, LUMINOSITY: {LUMINOSITY},
               'dimensionless': {NORMALIZED}}


def _select_models(index, parameters, teff_column, n_models):
    """Positions of models spread over teff, each near the median of the
    other parameters."""
    teffs = np.unique(index[teff_column].astype(float))
    picks = teffs[np.unique(np.round(
        np.linspace(0, len(teffs) - 1, min(n_models, len(teffs)))).astype(int))]
    others = [p for p in parameters if p != teff_column and p in index]
    positions = []
    for teff in picks:
        rows = np.where(np.isclose(index[teff_column].astype(float), teff))[0]
        distance = np.zeros(len(rows))
        for name in others:
            values = index[name].astype(float).values
            spread = np.ptp(values) or 1.0
            distance += ((values[rows] - np.median(values)) / spread) ** 2
        positions.append(int(rows[np.argmin(distance)]))
    return positions


def check_flux_scale(path, n_models=3, teff_column='teff'):
    """Compare the fluxes of a grid with the surface flux of a blackbody.

    Reads ``n_models`` model spectra, spread from the lowest to the highest
    effective temperature of the grid, without loading the whole flux
    array.

    Parameters
    ----------
    path : str or pathlib.Path
        Grid file.
    n_models : int
        Number of models to check.
    teff_column : str
        Index column holding the effective temperature in K.

    Returns
    -------
    FluxScaleCheck
    """
    path = Path(path).expanduser()
    meta = _read_grid_meta(path)
    index = pd.read_hdf(path, 'index')
    if teff_column not in index:
        raise ValueError('The grid index has no {0!r} column; give the '
                         'effective temperature column as teff_column'
                         .format(teff_column))
    parameters = list(meta['parameters'])
    unit = u.Unit(meta['flux_unit'])
    label_kind = _label_kind(unit)
    wavelength = u.Quantity(pd.read_hdf(path, 'wavelength').values[:, 0],
                            meta['wavelength_unit']).to_value(u.AA)

    models = []
    with h5py.File(path, 'r') as handle:
        for position in _select_models(index, parameters, teff_column,
                                       n_models):
            flux = np.asarray(handle['fluxes'][position], dtype=float)
            flux = _to_label_cgs(flux, unit, label_kind)
            teff = float(index[teff_column].iloc[position])
            finite = np.isfinite(flux)
            if finite.sum() < 2:
                models.append(ModelFluxCheck(
                    parameters={name: float(index[name].iloc[position])
                                for name in parameters if name in index},
                    ratio=float('nan'), median_flux=float('nan'),
                    interpretation=INVALID, implied_radius=None,
                    n_nonfinite=int((~finite).sum()), min_flux=float('nan')))
                continue
            ratio = (trapezoid(flux[finite], wavelength[finite])
                     / trapezoid(blackbody_surface_flux(wavelength[finite], teff),
                                 wavelength[finite]))
            median_flux = float(np.median(flux[finite]))
            min_flux = float(np.min(flux[finite]))
            n_nonfinite = int((~finite).sum())
            if n_nonfinite or min_flux < 0:
                interpretation = INVALID
            else:
                interpretation = _interpret(ratio, median_flux,
                                            float(np.max(flux)))
            radius = None
            if interpretation == LUMINOSITY:
                radius = float(np.sqrt(ratio / (4 * np.pi))
                               / const.R_sun.cgs.value)
            models.append(ModelFluxCheck(
                parameters={name: float(index[name].iloc[position])
                            for name in parameters if name in index},
                ratio=float(ratio), median_flux=median_flux,
                interpretation=interpretation, implied_radius=radius,
                n_nonfinite=n_nonfinite, min_flux=min_flux))

    kinds = sorted({model.interpretation for model in models
                    if model.interpretation != INVALID})
    if not kinds:
        interpretation = 'no valid model among those checked'
    elif len(kinds) == 1:
        interpretation = kinds[0]
    else:
        interpretation = 'mixed: ' + ', '.join(kinds)
    if len(kinds) != 1 or kinds[0] == UNRECOGNIZED or label_kind == 'other':
        consistent = None
    else:
        consistent = kinds[0] in _CONSISTENT.get(label_kind, set())
    return FluxScaleCheck(
        path=str(path), flux_unit=str(meta['flux_unit']),
        label_kind=label_kind, wavelength_min=float(wavelength[0]),
        wavelength_max=float(wavelength[-1]), models=models,
        interpretation=interpretation, consistent_with_label=consistent)


def format_flux_scale_check(check):
    """Text report of a ``FluxScaleCheck``."""
    lines = ['Grid: {0}'.format(check.path),
             'Flux unit label: {0} ({1})'.format(check.flux_unit,
                                                 check.label_kind),
             'Wavelength range: {0:.1f} to {1:.1f} Angstrom'.format(
                 check.wavelength_min, check.wavelength_max),
             'Grid flux / blackbody surface flux, integrated over the '
             'wavelength range:']
    for model in check.models:
        parameters = ', '.join('{0}={1:g}'.format(name, value)
                               for name, value in model.parameters.items())
        line = '  {0}: ratio {1:.3g}, {2}'.format(parameters, model.ratio,
                                                   model.interpretation)
        if model.implied_radius is not None:
            line += ', implied radius {0:.3g} R_sun'.format(
                model.implied_radius)
        if model.interpretation == INVALID:
            line += ' ({0} non-finite pixels, minimum flux {1:.3g})'.format(
                model.n_nonfinite, model.min_flux)
        lines.append(line)
    lines.append('Fluxes: {0}'.format(check.interpretation))
    n_invalid = sum(model.interpretation == INVALID for model in check.models)
    if n_invalid:
        lines.append('{0} of the {1} models checked are invalid and were not '
                     'used to classify the fluxes.'.format(n_invalid,
                                                           len(check.models)))
    if check.consistent_with_label is True:
        lines.append('The fluxes agree with the flux unit label.')
    elif check.consistent_with_label is False:
        lines.append('The fluxes do not agree with the flux unit label. '
                     'Operations that need absolute fluxes, such as '
                     'Distance and the magnitudes from Photometry, give '
                     'wrong results with this grid.')
    else:
        lines.append('The agreement with the flux unit label could not be '
                     'decided.')
    return '\n'.join(lines)
