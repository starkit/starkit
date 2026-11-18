********************
Synthetic Photometry
********************

``starkit.base.operations.imager.Photometry`` converts a spectrum to
magnitudes in a set of filters. The filters are curves from the SVO Filter
Profile Service, curves read from text files, or curves given as arrays. The
magnitudes are in the Vega or the AB system. The flux is taken to be
:math:`f_\lambda` in erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1` at the
observer; read :doc:`io/flux_units` before comparing the magnitudes of a grid
with observed ones.

No package beyond the StarKit requirements is needed; the Vega spectrum for
Vega magnitudes is included with StarKit.

Quick start
===========

Filters are named by their SVO filter IDs, as listed at
http://svo2.cab.inta-csic.es/theory/fps/::

    from starkit.gridkit import load_grid
    from starkit.base.operations.imager import Photometry

    grid = load_grid('phoenix_t4000_10000_w3000_9000_r3000.h5')
    model = grid | Photometry.from_grid(
        grid, filter_set=['SLOAN/SDSS.g', 'SLOAN/SDSS.r', 'Generic/Johnson.V'])
    model.teff_0, model.logg_0, model.mh_0 = 6000, 4.5, 0.0
    magnitudes = model()

``assemble_model`` adds the same operation when given ``filter_set`` and,
optionally, ``mag_type``::

    from starkit import assemble_model

    model = assemble_model(grid, filter_set=['SLOAN/SDSS.g', 'SLOAN/SDSS.r'],
                           mag_type='ab', vrad=0.)

Filter curves
=============

``filter_set`` is a list whose items are SVO filter IDs, ``Filter`` objects,
or both. ``starkit.filters.Filter`` holds a curve made in one of three ways::

    from astropy import units as u
    from starkit.filters import Filter

    # From the SVO Filter Profile Service
    ks = Filter.from_svo('2MASS/2MASS.Ks')

    # From a text file with two columns, wavelength and transmission;
    # lines starting with # are skipped
    nb217 = Filter.read('nb217.txt', wavelength_unit=u.um)

    # From arrays; wavelengths without a unit are taken as Angstrom
    box = Filter([4000., 4010., 5990., 6000.], [0., 1., 1., 0.], name='box')

    photometry = Photometry(filter_set=[ks, nb217, 'SLOAN/SDSS.r'])

Only the shape of a curve matters, so a transmission in percent gives the same
magnitudes as one in fractions. The wavelengths must increase. Between its
points a curve is interpolated linearly, and it is zero outside them.

``Filter.from_svo`` downloads the curve on first use and keeps it in the
astropy download cache, so later uses need no network connection. For
example, SVO has 11 Keck/NIRC2 curves, taken from the Keck NIRC2 filter
page: ``Keck/NIRC2.J``, ``.H``, ``.K``, ``.Kp``, ``.Ks``, ``.Lp`` and ``.Ms``,
and the narrow-band ``.Hcont``, ``.Fell`` (Fe II), ``.Brgamma`` and
``.Kcont``. Other NIRC2 filters can be read with ``Filter.read`` from a
tabulated curve. IDs in the
form used by wsynphot, such as ``'2MASS/2MASS/Ks'``, are also accepted.
wsynphot ``FilterSet`` and ``FilterCurve`` objects are accepted in
``filter_set``, without wsynphot being installed.

Detector type
=============

A response curve describes either the relative number of photons detected
at each wavelength or the relative energy detected. The mean flux density in
the filter is

.. math::

    \langle f_\lambda \rangle = \frac{\int f_\lambda S(\lambda) \lambda\,d\lambda}
                                     {\int S(\lambda) \lambda\,d\lambda}

for a photon-counting response :math:`S`, and

.. math::

    \langle f_\lambda \rangle = \frac{\int f_\lambda S'(\lambda)\,d\lambda}
                                     {\int S'(\lambda)\,d\lambda}

for an energy response :math:`S'`. The two describe the same system when
:math:`S' = S \lambda` (Bessell & Murphy 2012, PASP 124, 140, Appendix A).
CCDs and infrared arrays count photons, so their response curves are
photon-counting. A filter transmission is the fraction of the light
transmitted at each wavelength, the same for photons and for energy, so a
measured transmission used with such a detector is also a photon-counting
response. For a photomultiplier read by integrating its current, the system
response is the filter transmission multiplied by the radiant response of the
photocathode, in mA/W, which is an energy response; Bessell & Murphy (2012,
Section 5.4) note that the response curves generally published before CCDs
were of this kind.

The ``detector`` parameter of the ``Filter`` constructor and of
``Filter.read`` is ``'photon'`` by default and ``'energy'`` for an energy
response. ``Filter.from_svo`` takes it from the ``DetectorType`` that SVO
gives for each curve, unless ``detector`` is given. SVO lists the generic
Johnson U, B, V, R and I curves (``Generic/Johnson.*``) and the GCPD V curves
measured with photomultipliers as energy responses, and the V curves of CCD
instruments that include the detector response, such as
``Paranal/OmegaCAM.V``, as photon-counting. For curves of the filter alone
the SVO labels are not consistent: ``OSN/Circ.Johnson_V`` is listed as an
energy response and ``OSN/Johnson.V`` as photon-counting. A measured filter
transmission used with a CCD or an infrared array is photon-counting. Using
the generic Johnson V curve with the wrong integral changes the magnitude of
a 6000 K blackbody by 0.011 mag.

Curves that are all called Johnson V differ by more than this. For a 3500 K
blackbody, the curves ``CAHA/Astralux.V_noatm`` and
``GCPD/Johnson.V_Landolt`` give magnitudes 0.103 and 0.089 mag fainter than
``Generic/Johnson.V``; at 6000 K the differences are 0.043 and 0.036 mag. To
compare with observed magnitudes, use the curve of the system the
observations are calibrated to.

Zero points
===========

``mag_type='vega'`` gives magnitudes relative to the spectrum of Vega, with
Vega at 0 mag in every filter. The spectrum is the CALSPEC spectrum
``alpha_lyr_stis_011.fits`` (Bohlin, Gordon & Tremblay 2014, PASP 126, 711),
from 900 Å to 300 µm, which is included with StarKit; it is the default Vega
spectrum of synphot 1.9.0. Its mean flux density is computed with the same
integral as for the model spectrum. StarKit keeps the spectrum for the rest
of the Python session; to use another CALSPEC spectrum, set
``starkit.filters.VEGA_FILE`` to its path or URL and call
``starkit.filters.vega_spectrum.cache_clear()``.

``mag_type='ab'`` gives :math:`m = -2.5 \log_{10} (\langle f_\nu \rangle /
3631\ {\rm Jy})`, with :math:`\langle f_\nu \rangle = \langle f_\lambda
\rangle \lambda_p^2 / c` and the pivot wavelength

.. math::

    \lambda_p^2 = \frac{\int S \lambda\,d\lambda}{\int S / \lambda\,d\lambda}
                = \frac{\int S'\,d\lambda}{\int S' / \lambda^2\,d\lambda}

(Bessell & Murphy 2012, equation A16). ``Filter.zero_point(mag_type)`` returns
the zero point of a filter in erg s\ :sup:`-1` cm\ :sup:`-2` Å\ :sup:`-1`.

Wavelength coverage
===================

The spectrum must cover the filter. If more than 1% of the response of a
filter, :math:`\int S \lambda\,d\lambda` or :math:`\int S'\,d\lambda`, lies
outside the wavelengths of the spectrum, ``Photometry`` stops with a
``ValueError`` that names the filter. ``Photometry.from_grid`` checks this
when it is created, against the wavelengths of the grid.

``Photometry.from_grid`` interpolates the curves onto the grid wavelengths
once and uses them for every spectrum. A wavelength shift applied by an
earlier operation, such as ``DopplerShift``, therefore does not move the
spectrum relative to the filters. ``Photometry`` created without grid
wavelengths interpolates the curves onto the wavelengths it receives at each
evaluation.

Observed magnitudes in a fit
============================

``PhotometryColorLikelihood`` compares the colors between consecutive
filters, which do not depend on the flux scale of the grid.
``RelativePhotometryLikelihood`` compares the magnitudes themselves. Both
take the observed magnitudes and their uncertainties, in the order of the
filters::

    from starkit.fitkit.likelihoods import PhotometryColorLikelihood

    likelihood = model | PhotometryColorLikelihood(
        [17.64, 17.22, 17.03], [0.035, 0.053, 0.066])

An object with the attributes ``magnitudes`` and ``magnitude_uncertainties``,
such as a wsynphot ``MagnitudeSet``, is also accepted.

Accuracy
========

The magnitudes from ``Photometry`` were compared with those of
``Observation.effstim`` in synphot 1.9.0, using the same SVO curves and Vega
spectrum, for two
spectra: the spectrum of the example grid
``phoenix_t4000_10000_w3000_9000_r3000.h5`` at 6000 K, log g = 4.5 and
[M/H] = 0, and a 3000 K blackbody tabulated from 0.9 to 6 µm.

==========================================  ========  ===========  ===========
Filters                                     Detector  Vega (mag)   AB (mag)
==========================================  ========  ===========  ===========
SLOAN/SDSS.g, SLOAN/SDSS.r (6000 K grid)    photon    < 0.000001   0.000071
Paranal/OmegaCAM.V (6000 K grid)            photon    < 0.000001   0.000071
Generic/Johnson.B, .V (6000 K grid)         energy    0.000005     0.000066
11 Keck/NIRC2 curves (3000 K blackbody)     photon    0.000035     0.000066
==========================================  ========  ===========  ===========

The columns give the largest absolute difference between StarKit and synphot.
Of the AB differences, 0.000066 mag is the difference between the AB zero
point of 3631 Jy used by StarKit and the 3630.78 Jy (AB magnitude defined by
the constant -48.60) used by synphot through astropy. StarKit's Vega zero points
for the NIRC2 curves, converted to :math:`f_\nu` at the pivot wavelength,
agree with the zero points listed by SVO to 0.069%.

Two settings were needed on the synphot side. synphot holds the end values of
a curve that does not end at zero constant at all wavelengths beyond it;
StarKit sets the curve to zero there. The curves were therefore given to
synphot with ``fill_value=0``; with synphot's default, its magnitudes of the
3000 K blackbody in the narrow-band NIRC2 filters change by up to 0.645 mag.
synphot also treats every curve as photon-counting, so the energy responses
were interpolated onto a fine sampling and then divided by the wavelength.
Dividing the Johnson curves, which are sampled every 100 to 200 Å, before
interpolating them changes the synphot magnitudes by up to 0.0013 mag.

Changes from StarKit 0.8 and wsynphot
=====================================

StarKit 0.8 and earlier computed the magnitudes with wsynphot. The changes
that affect results are:

- The Vega spectrum. wsynphot used the CALSPEC model ``alpha_lyr_mod_002.fits``.
  Vega magnitudes computed with it differ from the present ones by
  -0.010 mag in ``Generic/Johnson.V``, -0.004 mag in ``SLOAN/SDSS.r`` and
  +0.020 mag in ``2MASS/2MASS.Ks`` (wsynphot minus present). Setting
  ``starkit.filters.VEGA_FILE`` to
  https://archive.stsci.edu/hlsps/reference-atlases/cdbs/calspec/alpha_lyr_mod_002.fits
  and calling ``starkit.filters.vega_spectrum.cache_clear()`` reproduces the
  wsynphot zero points.
- The detector type. wsynphot from October 2017 to January 2022 used the
  photon-counting integral for every curve, including the curves SVO lists
  as energy responses, such as the generic Johnson filters; earlier versions
  used the energy integral for every curve. Since February 2022 it used the
  SVO detector type, as StarKit does now.
- Observed magnitudes. The likelihoods take arrays of magnitudes and
  uncertainties instead of a wsynphot ``MagnitudeSet``, which is still
  accepted.
