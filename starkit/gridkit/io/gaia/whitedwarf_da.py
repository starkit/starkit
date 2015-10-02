import numpy as np
from scipy.interpolate import interp1d

from starkit.gridkit.io.gaia.base import GaiaSpecLibIO
from starkit.gridkit.io.process import BaseProcessGrid

class WhiteDwarfDAProcess(BaseProcessGrid):
    pass




class WhiteDwarfDAIO(GaiaSpecLibIO):

    R=10000

    def __init__(self, fname):
        super(WhiteDwarfDAIO, self).__init__(fname)
        self.process = WhiteDwarfDAProcess(self.get_grid_wavelength(), self.R, sampling=4)

    def get_wavelength(self):
        return self.process.output_wavelength * self.wavelength_unit

    def get_fluxes(self):
        fluxes = super(WhiteDwarfDAIO, self).get_fluxes()
        proc_fluxes = []

        for flux in fluxes:
            proc_fluxes.append(self.process.process(flux))
        return np.array(proc_fluxes) * self.flux_unit


    def to_hdf(self, fname, index_columns=['teff', 'logg'], clobber=False):
        super(WhiteDwarfDAIO, self).to_hdf(fname, index_columns,
                                           clobber=clobber)



