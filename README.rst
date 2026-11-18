StarKit
=======

Installation
************

StarKit requires Python 3.10 or later. The simplest setup is the conda
environment in ``environment.yml``. It uses Python 3.10 and contains StarKit
(installed in editable mode), its optional features, the MultiNest library
and pytest::

   git clone https://github.com/starkit/starkit.git
   cd starkit
   conda env create -f environment.yml
   conda activate starkit
   pytest

To install with pip instead::

   python -m pip install -e '.[all,test]'
   pytest

Optional extras currently available::

   all         grid_io, multinest, parallel and photometry together
   grid_io     Grid preparation helpers (progressbar2, sqlalchemy)
   multinest   MultiNest fitting and corner plots (pymultinest, corner)
   parallel    Distributed fitting helpers (ipyparallel, psutil)
   photometry  Vega magnitudes in synthetic photometry (synphot)
   test        pytest
   docs        Sphinx, its extensions and IPython, for building the documentation

``pymultinest`` is a Python interface to the native MultiNest library. The
conda environment includes the library. With pip, the library must be built
separately, and its directory must be on ``LD_LIBRARY_PATH`` (Linux) or
``DYLD_LIBRARY_PATH`` (macOS).

If you plan to use OpenMPI, also install ``mpi4py``::

   python -m pip install mpi4py

Synthetic photometry (``Photometry``) downloads filter curves from the SVO
Filter Profile Service. Vega magnitudes need the ``photometry`` extra
(synphot); AB magnitudes need no extra package.


Grid flux units
***************

Check the units of a grid before using it for absolute fluxes, such as with
``Distance``, the magnitudes from ``Photometry`` or
``RelativePhotometryLikelihood``. The ``flux_unit`` label of many grids is
wrong: StarKit labeled every grid it wrote between December 2017 and January
2025 ``erg/s/angstrom`` (luminosity density), although many of these grids
hold surface flux, in erg/(s cm2 Angstrom). ``Distance`` is correct only for
luminosity density; for surface flux, the flux at the observer is the surface
flux times (R/d)^2, so absolute fluxes also need the radius R of the star.
To check a grid::

   from starkit.gridkit.flux_check import check_flux_scale, format_flux_scale_check

   print(format_flux_scale_check(check_flux_scale('grid.h5')))

Spectral fits that normalize the model to the observed spectrum, and fits to
colors, do not depend on the flux scale. ``BUGS.md`` and the documentation
page on grid flux units give the details and a workaround.


Example publications that use StarKit
**************************************

- Do, Tuan; Kerzendorf, Wolfgang; Konopacky, Quinn; Marcinik, Joseph M.; Ghez, Andrea; Lu, Jessica R.; Morris, Mark R., `Super-solar Metallicity Stars in the Galactic Center Nuclear Star Cluster: Unusual Sc, V, and Y Abundances <https://ui.adsabs.harvard.edu/#abs/2018ApJ...855L...5D/abstract>`_
- Feldmeier-Krause, A.; Kerzendorf, W.; Neumayer, N.; Schödel, R.; Nogueras-Lara, F.; Do, T.; de Zeeuw, P. T.; Kuntschner, H., `KMOS view of the Galactic Centre - II. Metallicity distribution of late-type stars <https://ui.adsabs.harvard.edu/#abs/2017MNRAS.464..194F/abstract>`_

A few test grids can be found at https://starkit.github.io/starkit/io/available_grids.html. If you use any of these grids, please make sure that the grid from which it is created (like Phoenix grid) is also cited.
