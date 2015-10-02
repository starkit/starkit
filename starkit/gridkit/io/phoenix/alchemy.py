from astropy.io import fits

from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm import backref, relationship
from sqlalchemy import Column, ForeignKey
from sqlalchemy import Integer, Float, String

from astropy import units as u
from astropy.io import fits
from specutils import Spectrum1D

from starkit.gridkit.io.phoenix import parameters

PhoenixBase = declarative_base()



import os
import numpy as np

class Spectrum(PhoenixBase):
    __tablename__ = 'spectra'

    id = Column(Integer, primary_key=True)
    fpath = Column(String)
    fname = Column(String)

    base_dir = None
    wavelength = None
    flux_unit = 'erg/s/angstrom'

    @property
    def full_path(self):
        if self.base_dir is None:
            raise AttributeError('base_dir attribute needs to be set')
        else:
            return os.path.join(self.base_dir, self.fpath, self.fname)

    def __repr__(self):
        return '<Spectrum {0}>'.format(self.fname)

    def get_spectrum1d(self):
        flux = self._read_flux()
        return Spectrum1D.from_array(self.wavelength, flux,
                                     unit=u.Unit(self.flux_unit))

    def _read_flux(self):
        return fits.getdata(self.full_path) * 1e-8 * (
            4 * np.pi * fits.getval(self.full_path, 'PHXREFF')**2)


class ParameterSetMixin(object):
    ___tablename__ = 'parameter_sets'

    @declared_attr
    def spectrum_id(cls):
        return Column(Integer, ForeignKey('spectra.id'))

    @declared_attr
    def spectrum(cls):
        return relationship(Spectrum, uselist=False,
                            backref=backref('parameter_set', uselist=False))

    @classmethod
    def from_file(cls, fname):
        param_dict = {}
        for param in cls.parameters:
            param_dict[param.name] = param.ingest(fname)
        return cls(**param_dict)

    def __repr__(self):
        param_strs = ['{0} = {1}'.format(param.name, getattr(self, param.name))
                      for param in self.parameters]
        return '<ParameterSet\n{0}\n>'.format('\n'.join(param_strs))


def make_parameter_set():
    type_dict = {'float':Float}
    parameter_classes = parameters.BasePhoenixParameter.__subclasses__()
    class_dict = {'parameters': [item() for item in parameter_classes]}
    class_dict['id'] = Column(Integer, primary_key=True)
    for param in parameter_classes:
        class_dict[param.name] = Column(type_dict[param.type])

    class_dict['__tablename__'] = 'parameter_sets'

    return type('ParameterSet', (ParameterSetMixin, PhoenixBase), class_dict)

ParameterSet = make_parameter_set()