*****************************
Working with the Phoenix Grid
*****************************

The Phoenix Grid is available via ftp at
`ftp://phoenix.astro.physik.uni-goettingen.de/HiResFITS/`_. The current
implementation of StarKit only works with the HiResFITS version of the
Phoenix grid. Download the entire grid or portions of it. Then generate a
database (it only saves the location of the files and the stellar parameters -
not the grid itself)::

    from starkit.gridkit.io.phoenix.base import PhoenixSpectralGridIO
    grid = PhoenixSpectralGridIO('sqlite:///phoenix.db3', base_dir='/media/data1/grids/phoenix')

where `phoenix.db3` is the filename for your database and `base_dir` needs to have
 all the `Z*` directories of the phoenix grid.

To load the phoenix grid into a database (it will only save the parameters and
a path to the file in the database)::

    grid.ingest()

This will take some time so get a cup of coffee.

After this the grid can be cut up and be made useable for working with Starkit
by writing it to HDF5::

    from starkit.gridkit.io.phoenix.base import PhoenixSpectralGridIO
    from starkit.gridkit.io.phoenix.alchemy import ParameterSet

    from astropy import units as u

    pgrid = PhoenixSpectralGridIO('sqlite:///phoenix.db3', '/media/data1/grids/phoenix')
    fluxes = pgrid.to_hdf('phoenix_mcsnr.h5', (ParameterSet.mh>-1.5, ParameterSet.teff.between(5000, 9000)),
    wavelength_range=(2000, 9000),R=10000, clobber=True)

The syntax for the parameter filter uses the sqlalchemy backend.