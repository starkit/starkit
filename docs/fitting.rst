*******
Fitting
*******

StarKit fits spectra with MultiNest, through the ``MultiNest`` class in
``starkit.fitkit.multinest``. MultiNest fitting needs the pymultinest package
and the compiled MultiNest library, and a default installation of StarKit
includes neither. The conda environment in ``environment.yml`` includes both.
With pip, install the ``multinest`` extra and build the library separately.
:doc:`installation` describes both methods.

MultiNest
---------

When using MultiNest please make sure to cite their works (arXiv:0809.3437,
arXiv:0704.3704, arXiv:1306.2144) as well as the wrapper PyMultiNest
( `Buchner et al. 2014, A&A <http://www.aanda.org/articles/aa/abs/2014/04/aa22971-13/aa22971-13.html>`_)
