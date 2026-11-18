# Known bugs

All entries were found on 2026-11-17 in StarKit 0.7 and are open.

1. [Grid flux units and operations that need absolute fluxes](#grid-flux-units-and-operations-that-need-absolute-fluxes)
2. [PHOENIX grids have half the resolving power in their metadata](#phoenix-grids-have-half-the-resolving-power-in-their-metadata)
3. [Pollux grid builder](#pollux-grid-builder)
4. [The loglikelihood column of MultiNest results holds -2 ln L](#the-loglikelihood-column-of-multinest-results-holds--2-ln-l)
5. [InstrumentDeltaLambdaConstant with numpy 2.5](#instrumentdeltalambdaconstant-with-numpy-25)
6. [Normalize warns about well-conditioned fits](#normalize-warns-about-well-conditioned-fits)
7. [Gaia and Munari grid builders](#gaia-and-munari-grid-builders)
8. [Bell characters in the PHOENIX and BOSZ bibtex](#bell-characters-in-the-phoenix-and-bosz-bibtex)

[Checked and correct](#checked-and-correct) and [Not yet checked](#not-yet-checked)
list related findings that are not bugs.

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

## PHOENIX grids have half the resolving power in their metadata

`PhoenixProcessGrid.interp_wavelength` first interpolates each spectrum onto
wavelengths with `pre_sampling` pixels per resolution element of the PHOENIX
spectra, then convolves with a Gaussian whose width in pixels it computes with
`R_sampling` instead of `pre_sampling`. The Gaussian is therefore
`R_sampling / pre_sampling` times the intended width: twice with the defaults,
`R_sampling=4` and `pre_sampling=2`. The metadata records the requested `R`.

Measured on a synthetic line at 5000 Angstrom with an intrinsic FWHM of
0.01 Angstrom, processed with `R=10000` and `R_sampling=4`: a Gaussian fit
gives a FWHM of 1.000 Angstrom (R = 5000) with `pre_sampling=2` and
0.500 Angstrom (R = 10000) with `pre_sampling=4`. The code dates from
2017-12-19 (f9a8ea8), so PHOENIX grids built with the defaults since then are
probably affected; their metadata does not record `pre_sampling`.

Fitted parameters that set the line widths (`vrot`, an instrumental `R`)
absorb the extra width, so they are biased; `InstrumentRConstant` computes its
kernel from the `R` in the metadata and leaves the model broader than
requested.

Workaround: build grids with `pre_sampling` equal to `R_sampling`.
Fix: compute the kernel width with `pre_sampling`, or better, from the
sampling of the interpolated wavelengths.

## Pollux grid builder

`POLLUXCMFGENProcessGrid` (`starkit/gridkit/io/pollux/process.py`) has three
problems.

- **Kernel width.** It assumes input spectra with R = 50000 sampled at
  2 pixels per resolution element, 3.0 km/s per pixel. The Pollux CMFGEN
  infrared spectra have 3.75 km/s per pixel (0.150 Angstrom at
  12000 Angstrom, 0.314 Angstrom at 25000 Angstrom), so the Gaussian is
  1.25 times too wide. A line with an intrinsic FWHM of 6 km/s processed with
  `R=5000` comes out with R = 4002. Fix: compute the input sampling from the
  wavelengths, as for PHOENIX.
- **Flux cache reader.** It imports `convert_bz2_memmap` from the BOSZ module,
  which builds the cache name by replacing `.bz2` and does nothing for a
  `.spec` file. Without existing caches, `load_flux` stops with
  `FileNotFoundError`. The caches must first be written with
  `starkit.gridkit.io.pollux.base.cache_pollux_grid`. Fix: import the
  Pollux reader.
- **Flux and wavelength lengths.** A Pollux grid built in September 2025 was built from
  caches holding the flux column without its first row (58638 values for
  58639 wavelengths), written by code that is not in the repository. Each
  flux was then paired with the next wavelength, which shifts the models by
  one input pixel, 3.75 km/s to the blue (3.3 km/s measured by
  cross-correlation), and fitted radial velocities by about +3.75 km/s.
  Nothing checks that a flux array and the wavelength array have the same
  length. Fix: check the lengths in `BaseProcessGrid.get_fluxes`, and rename
  the caches (`.v2.npy`) so that caches from older code are not reused.

## The loglikelihood column of MultiNest results holds -2 ln L

`MultiNestResult.read_posterior_data` names the second column of the MultiNest
`*_.txt` file `loglikelihood`. MultiNest writes -2 ln L there: for one fit,
the row with the lowest value had 570.66 in that column and ln L = -285.33.
`MultiNestResult.maximum` takes the minimum of the column and so returns the
sample with the highest likelihood, which is correct. Fix: rename the column,
or store -0.5 times the value.

## InstrumentDeltaLambdaConstant with numpy 2.5

`InstrumentDeltaLambdaConstant.evaluate` (the `delta_lambda` keyword of
`assemble_model`) passes its one-element parameter array to `np.arange`, and
with numpy 2.5 the evaluation stops with `TypeError: only 0-dimensional arrays
can be converted to Python scalars`. Fix: convert the parameter to a scalar,
as the stellar operations do with `_as_scalar_float`.

## Normalize warns about well-conditioned fits

`Normalize.evaluate` warns `The fit may be poorly conditioned` when the rank of
its least-squares problem differs from `npol`. A polynomial of degree `npol`
has `npol + 1` coefficients, so the warning appears for every full-rank fit,
and not for a fit of rank exactly `npol`. With `npol=3`, it was emitted for
ranks 4 and 2 and not for rank 3. Fix: compare the rank with `npol + 1`.

## Gaia and Munari grid builders

`GaiaSpecLibIO` and `Munari1A` write grids with
`BaseSpectralGridIO.to_hdf` (`starkit/gridkit/io/base.py`), which stores the
fluxes and wavelengths in the older layout without a `meta` table;
`load_grid` cannot read the result. `Munari1A.to_hdf` also passes a
`default_values` argument that `BaseSpectralGridIO.to_hdf` does not accept.
Importing `starkit.gridkit.io.base` requires sqlalchemy (the `grid_io`
extra). Fix: move these builders to `BaseProcessGrid`, or remove them.

## Bell characters in the PHOENIX and BOSZ bibtex

`phoenix_bibtex` and `bosz_bibtex` are ordinary strings containing `{\aap}` and
`{\aj}`, in which Python reads `\a` as the bell character (U+0007). The
`bibtex` metadata of PHOENIX and BOSZ grids therefore has a control character
where the journal macro should be. Fix: make them raw strings, as
`gotberg23_bibtex` is.

## Checked and correct

- **Pollux wavelength convention.** The `.spec.txt` header of
  the Pollux CMFGEN infrared spectra says `column1 = 'wavelength' /air
  wavelength`, but the line cores of H I Brackett 4-7 to 4-11 in five models
  (12020 to 24030 K) of a Pollux grid lie at their vacuum wavelengths
  (Rydberg formula), with a median offset of -8.4 km/s that includes the
  -3.75 km/s shift described above; air wavelengths would put them at
  -85.6 km/s. The builder's `wavelength_type='vacuum'` is correct and should
  not be taken from the header.

## Not yet checked

- **CMFGEN builder** (`starkit/gridkit/io/cmfgen`). Its cache reader treats
  the first line of each `.tsv` file as a header, so it drops a row if the
  files have none, and it declares vacuum wavelengths and erg/(s cm2
  Angstrom) without converting anything. Both depend on the format of the
  input files, which has not been checked.
- **BOSZ builder.** It assumes input spectra with R = 300000 sampled at one
  pixel per resolution element; the sampling of the BOSZ files has not been
  checked against this.
