***************
Getting Started
***************

.. _working-with-grid:

Working with the grid
^^^^^^^^^^^^^^^^^^^^^

After successful creation of an HDF5 Grid (please see :ref:`io` for details)
one can load the grid using::

    >>> from starkit.gridkit import load_grid
    >>> grid = load_grid('mygrid.h5')
    >>> grid
    <SpectralGrid(teff=6200.0, logg=4.0, mh=0.5)>

The values for `teff`, `logg` and `mh` are the defaults for the grid, but have
no other special meaning and can be changed.
StarKit heavily uses astropy.modeling please refer to its
`reference guide <http://astropy.readthedocs.org/en/latest/modeling/>`_.

One can then just evaluate the grid at a certain parameter point::

    >>> grid.teff = 5780.
    >>> grid.logg = 4.4
    >>> grid.mh = 0.0
    >>> wave, flux = grid()
    >>> plot(wave, flux)

StarKit works with a so called plugin-architecture that can modify the generated
spectrum. For example, one can use the doppler shift plugin::

    >>> from starkit.base.operations.stellar import DopplerShift, RotationalBroadening
    >>> doppler = DopplerShift(vrad=0)
    >>> doppler.vrad = 20
    >>> new_wave, new_flux = doppler(wave, flux)

These plugins can then be concatenated (again see
`reference guide <http://astropy.readthedocs.org/en/latest/modeling/>`_)::

    >>> my_model = grid | doppler
    >>> my_model
    <CompoundModel0(teff_0=5780.0, logg_0=4.4, mh_0=0.0, vrad_1=20.0)>

You can then add further plugins to the model::

    >>> from starkit.base.operations.stellar import RotationalBroadening, CCM89Extinction
    >>> rotation = RotationalBroadening.from_grid(grid, vrot=200) # it is imperative to load from_grid here as it needs to
    >>> extinction = CCM89Extinction(a_v=1.0, r_v=3.1)
    >>> my_model = grid | rotation | doppler | extinction
    >>> my_model
    <CompoundModel3(teff_0=2300.0, logg_0=0.0, mh_0=0.5, alpha_0=0.0, vrot_1=200.0, limb_darkening_1=0.6, vrad_2=200.0, a_v_3=1.0, r_v_3=3.1)>
    >>> my_model.param_names
    (u'teff_0',
     u'logg_0',
     u'mh_0',
     u'alpha_0',
     u'vrot_1',
     u'limb_darkening_1',
     u'vrad_2',
     u'a_v_3',
     u'r_v_3')
     >>> my_model.vrot_1 = 300

There is a convenience function to assemble models::

    >>> from starkit import assemble_model
    >>> my_model = assemble_model(grid, vrad=20, a_v=0.3, vrot=300)
    >>> my_model.__class__
    <class 'starkit.base.assemble_model.CompoundModel5'>
    Name: CompoundModel5
    Inputs: ()
    Outputs: ('wavelength', 'flux')
    Fittable parameters: (u'teff_0', u'logg_0', u'mh_0', u'alpha_0', u'vrot_1', u'limb_darkening_1', u'vrad_2', u'a_v_3', u'r_v_3')
    Expression: [0] | [1] | [2] | [3]
    Components:
        [0]: <SpectralGrid(teff=2300.0, logg=0.0, mh=0.5, alpha=0.0)>

        [1]: <RotationalBroadening(vrot=300.0, limb_darkening=0.6)>

        [2]: <DopplerShift(vrad=20.0)>

        [3]: <CCM89Extinction(a_v=0.3, r_v=3.1)>

Working with grid and spectra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The next step is to work with grid and spectra.