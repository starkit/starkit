from astropy.io import fits
from abc import abstractproperty, abstractmethod, ABCMeta



class BasePhoenixData(object):

    __metaclass__ = ABCMeta

    @abstractproperty
    def name(self):
        raise NotImplementedError

    @abstractmethod
    def ingest(self, fname):
        raise NotImplementedError('needs to be implemented by subclasses')

    def __repr__(self):
        return "<Parameter {0}>".format(self.name)

class BasePhoenixParameter(BasePhoenixData):
    __metaclass__ = ABCMeta

class BasePhoenixMetaData(BasePhoenixData):
    __metaclass__ = ABCMeta


class Teff(BasePhoenixParameter):
    name = 'teff'
    type = 'float'

    def ingest(self, fname):
        return fits.getval(fname, 'PHXTEFF')


class Logg(BasePhoenixParameter):
    name = 'logg'
    type = 'float'

    def ingest(self, fname):
        return fits.getval(fname, 'PHXLOGG')

class MH(BasePhoenixParameter):
    name = 'mh'
    type = 'float'

    def ingest(self, fname):
        return fits.getval(fname, 'PHXM_H')

class Alpha(BasePhoenixParameter):
    name = 'alpha'
    type = 'float'

    def ingest(self, fname):
        return fits.getval(fname, 'PHXALPHA')


