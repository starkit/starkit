# Known bugs

## Grid flux units and operations that need absolute fluxes

Status: open. Found 2026-11-17 in StarKit 0.7.

### Summary

The `flux_unit` label of many grids does not describe their fluxes, and
`Distance`, the only operation that scales a model to the observer, is correct
only for grids of luminosity density, while the grid builders write surface
flux. Any result that depends on absolute fluxes (`Distance`, magnitudes from
`Photometry`, `RelativePhotometryLikelihood`) is wrong for most grids, by
about 57 mag for a star of one solar radius. Spectral fits that normalize the
model to the observed spectrum (`normalize_npol`), and fits to colors
(`PhotometryColorLikelihood`), do not depend on the flux scale and are not
affected.

### Flux conventions

- Surface flux, F_lambda, in erg/(s cm2 Angstrom). The flux observed at
  distance d from a star of radius R is F_lambda (R/d)^2.
- Luminosity density, L_lambda = 4 pi R^2 F_lambda, in erg/(s Angstrom). The
  observed flux is L_lambda / (4 pi d^2).
- Continuum-normalized flux, dimensionless, with no absolute scale.

### Labels written by StarKit

| Code | Fluxes written | `flux_unit` label |
|---|---|---|
| `BaseProcessGrid.get_meta`, 2017-12-20 (ff2c6dc) to 2025-01-15 (8c5fc1f) | whatever the builder produced | `erg/s/angstrom` for every grid |
| `BaseProcessGrid.get_meta`, since 2025-01-15 (8c5fc1f) | whatever the builder produced | `erg/s/cm^2/angstrom` for every grid, including telluric transmission |
| `PhoenixProcessGrid` before 2018-10-19 (5566c62) | luminosity density, with the radius `PHXREFF` of each PHOENIX model | |
| `PhoenixProcessGrid` since 2018-10-19 | surface flux | |
| `BOSZProcessGrid` before 2018-05-29 (1e65adb) | surface flux divided by pi | |
| `BOSZProcessGrid` since 2018-05-29 | surface flux | |

The label is therefore wrong for PHOENIX grids built with code from
2018-10-19 to 2025-01-15 and for BOSZ grids built before 2025-01-15. The
fluxes written by the CMFGEN and Gotberg 2023 builders have not been checked;
a Pollux grid built in September 2025 holds surface flux, as labeled. The date
of a grid file does not show which code built it. Grids built outside StarKit
can carry any label; a grid of continuum-normalized spectra labeled
`erg/s/angstrom` has been seen.

The test grid `docs/examples/phoenix_t4000_10000_w3000_9000_r3000.h5` holds
surface flux and is labeled `erg/s/angstrom`: its spectra integrate to 0.83,
0.99 and 0.90 times the surface flux of a blackbody at 4000, 6200 and
10000 K over 3000 to 9000 Angstrom.

### Effect on the operations

- `Distance` divides the flux by 4 pi d^2.
  - On a grid labeled as surface flux it stops with `UnitConversionError`.
  - On a surface-flux grid labeled `erg/s/angstrom` it runs and returns a flux
    too low by 4 pi R^2 (R in cm), 57.0 mag for R = 1 R_sun. The saved output
    of the example notebook "Starkit Photometry and Spectroscopy Combined
    Model" (Johnson U, B and V of 57.27 to 57.35 mag for a 9940 K model at
    3.5 pc) had this error.
  - Before StarKit 0.7, `Distance.from_grid` read `grid.flux_unit`, which
    grids do not have, so it stopped with `AttributeError` unless that
    attribute was set by hand. Since 0.7 it reads `grid.fluxes.unit`, so the
    wrong result above is reached without any warning.
- `Photometry` computes magnitudes assuming erg/(s cm2 Angstrom) at the
  observer. For a surface-flux grid without scaling, they are magnitudes of
  the stellar surface.
- `RelativePhotometryLikelihood` compares magnitudes, so its likelihood is
  wrong by the same amount.

### Checking a grid

`starkit.gridkit.flux_check.check_flux_scale(path)` compares a few model
spectra with the surface flux of a blackbody at the same effective
temperature and reports whether the grid holds surface flux, surface flux
divided by pi, luminosity density (with the implied radius) or
continuum-normalized spectra, and whether that agrees with the label. See
`docs/io/flux_units.rst`.

### Workaround

For a surface-flux grid, scale the flux by (R/d)^2 instead of using
`Distance`:

```python
from astropy import units as u
from astropy.modeling.models import Identity, Scale

dilution = Scale(float((1.0 * u.R_sun / (10 * u.pc)).to(1) ** 2))
model = assemble_model(grid, vrad=0.) | (Identity(1) & dilution)
```

R and d cannot be fitted separately from the factor (R/d)^2.

### Planned fix

1. An operation for surface-flux grids with the parameters `radius` and
   `distance`, selected by `Distance.from_grid` and `assemble_model` from the
   grid's flux unit, and a clear error when `distance` is used on a
   surface-flux grid without a radius.
2. Each grid builder declares the unit its `load_flux` returns, and
   `get_meta` writes that unit instead of a fixed string; dimensionless grids
   (telluric transmission, normalized spectra) are labeled dimensionless.
3. A `flux_unit=` argument to `load_grid` to override a wrong label, and a
   warning when a grid with `format_version` 1.0 is labeled `erg/s/angstrom`.
4. `format_version` 1.1 for grids written with the corrected labels.
