*******************
GridKit File Format
*******************

All StarKit input grids are stored in an HDF5 file. It consists of three parts:

1) The index is stored as a pandas DataFrame (through HDFStore) on the key 'index'.


2) The flux array is stored as a 2D array using h5py under the key 'fluxes'


The metadata is stored as a pandas Series (through HDFStore) on the key 'meta'

The current version of StarKit requires the following metadata:

===================  ======================================
flux_unit            astropy compliant unit
grid_bibtex          bibtex for stellar grid
grid_name            Grid Name
resolution_profile   only 'gaussian' supported at this time
resolution_sampling  how many pixels per resolution element
resolution_type      r_const or delta_lambda_const
uuid                 UUID4
wavelength_unit      astropy compliant unit
===================  ======================================