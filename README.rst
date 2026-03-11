StarKit
=======

Installation
************

Modernized installation procedure in this fork, to support python 3.12 more cleanly::

   # for Linux and Python 3.12
   git clone https://github.com/thorsbro/starkit.git
   cd starkit
   git remote add upstream https://github.com/starkit/starkit.git
   python -m pip install -e .

Install the most common optional features as well::

   python -m pip install -e '.[all]'

Optional extras currently available::

   all         Common optional features (ipyparallel, jbopt, progressbar2, psutil)
   grid_io     Grid preparation helpers (progressbar2)
   parallel    Distributed fitting helpers (ipyparallel, psutil)
   optimizers  Additional optimizer support (jbopt)
   photometry  Photometry support (wsynphot)

If you plan to use OpenMPI, also install ``mpi4py``::

   python -m pip install mpi4py

The photometry extra depends on ``wsynphot``, which may need to be installed
separately depending on how that project is distributed.


Example publications that use StarKit
**************************************

- Do, Tuan; Kerzendorf, Wolfgang; Konopacky, Quinn; Marcinik, Joseph M.; Ghez, Andrea; Lu, Jessica R.; Morris, Mark R., `Super-solar Metallicity Stars in the Galactic Center Nuclear Star Cluster: Unusual Sc, V, and Y Abundances <https://ui.adsabs.harvard.edu/#abs/2018ApJ...855L...5D/abstract>`_
- Feldmeier-Krause, A.; Kerzendorf, W.; Neumayer, N.; Schödel, R.; Nogueras-Lara, F.; Do, T.; de Zeeuw, P. T.; Kuntschner, H., `KMOS view of the Galactic Centre - II. Metallicity distribution of late-type stars <https://ui.adsabs.harvard.edu/#abs/2017MNRAS.464..194F/abstract>`_

A few test grids can be found at https://starkit.github.io/starkit/io/available_grids.html. If you use any of these grids, please make sure that the grid from which it is created (like Phoenix grid) is also cited.
