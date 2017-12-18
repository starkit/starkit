.. _io:

****************
HDF5 description
****************

All StarKit input grids are stored in an HDF5 file. It consists of two portions: The index is stored
as a pandas DataFrame (through HDFStore) on the key 'index'. The flux array is stored as a 2D array
