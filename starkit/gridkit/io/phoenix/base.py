from glob import glob
from astropy.io import fits
import pandas as pd
import numpy as np
from progressbar import ProgressBar


phoenix_bibtex = """
@ARTICLE{2013A&A...553A...6H,
   author = {{Husser}, T.-O. and {Wende-von Berg}, S. and {Dreizler}, S. and 
	{Homeier}, D. and {Reiners}, A. and {Barman}, T. and {Hauschildt}, P.~H.
	},
    title = "{A new extensive library of PHOENIX stellar atmospheres and synthetic spectra}",
  journal = {\aap},
archivePrefix = "arXiv",
   eprint = {1303.5632},
 primaryClass = "astro-ph.SR",
 keywords = {stars: atmospheres, convection, stars: late-type},
     year = 2013,
    month = may,
   volume = 553,
      eid = {A6},
    pages = {A6},
      doi = {10.1051/0004-6361/201219058},
   adsurl = {http://adsabs.harvard.edu/abs/2013A%26A...553A...6H},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
"""

phoenix_meta = {'bibtex':phoenix_bibtex,
                'parameters':['teff', 'logg', 'mh', 'alpha'],
                'wavelength_unit':'Angstrom',
                'wavelength_type':'vacuum',
                'flux_unit': 'erg/s/cm^2/angstrom'}


def make_raw_index():
    """
    Read all Phoenix files and generate a raw index with filename association.

    Returns
    -------
        phoenix_index : pd.DataFrame
    """
    all_fnames = glob('PHOENIX-ACES-AGSS-COND-2011/Z*/*.fits')
    phoenix_index = pd.DataFrame(index=np.arange(len(all_fnames)), columns=['teff', 'logg', 'mh', 'alpha', 'filename'])
    print("Reading Phoenix grid...")
    progressbar = ProgressBar(max_value=len(all_fnames))
    for i, fname in progressbar(enumerate(all_fnames)):
        spec_header = fits.getheader(fname)
        phoenix_index.iloc[i] = (spec_header['PHXTEFF'], spec_header['PHXLOGG'], spec_header['PHXM_H'],
                                 spec_header['PHXALPHA'], fname)
    return phoenix_index


def make_grid_info(fname):
    """
    Make the HDF5 Grid Info file

    Parameters
    ----------
    fname: str

    """

    raw_index = make_raw_index()
    wavelength = fits.getdata('WAVE_PHOENIX-ACES-AGSS-COND-2011.fits')

    with pd.HDFStore(fname) as fh:
        fh['index'] = raw_index
        fh['wavelength'] = pd.DataFrame(wavelength)
        fh['meta'] = pd.Series(phoenix_meta)




