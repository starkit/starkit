import re
import os

from glob import glob

from astropy.io import fits
import pandas as pd
import numpy as np
from progressbar import ProgressBar



bosz_bibtex = """
@ARTICLE{2017AJ....153..234B,
   author = {{Bohlin}, R.~C. and {M{\'e}sz{\'a}ros}, S. and {Fleming}, S.~W. and 
	{Gordon}, K.~D. and {Koekemoer}, A.~M. and {Kov{\'a}cs}, J.},
    title = "{A New Stellar Atmosphere Grid and Comparisons with HST/STIS CALSPEC Flux Distributions}",
  journal = {\aj},
archivePrefix = "arXiv",
   eprint = {1704.00653},
 primaryClass = "astro-ph.SR",
 keywords = {stars: atmospheres, stars: fundamental parameters, techniques: spectroscopic},
     year = 2017,
    month = may,
   volume = 153,
      eid = {234},
    pages = {234},
      doi = {10.3847/1538-3881/aa6ba9},
   adsurl = {http://adsabs.harvard.edu/abs/2017AJ....153..234B},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

"""

bosz_meta = {'bibtex':bosz_bibtex,
                'parameters':['teff', 'logg', 'mh', 'ch', 'alpha', 'rot', 'micro'],
                'wavelength_unit':'Angstrom',
                'wavelength_type':'vacuum'}


def make_raw_index():
    """
    Read all Phoenix files and generate a raw index with filename association.

    Returns
    -------
        bosz_index : pd.DataFrame
    """
    all_fnames = glob('ascii/insbroad_300000/*/*/*/*.asc.bz2')
    bosz_index = pd.DataFrame(index=np.arange(len(all_fnames)), columns=['teff', 'logg', 'mh', 'alpha', 'ch',
                                                                            'rot', 'micro', 'res', 'filename'])
    print "Reading Phoenix grid..."
    bar = ProgressBar(max_value=len(all_fnames))
    pattern = re.compile('am(p|m)(\d+)c(p|m)(\d+)+o(p|m)(\d+)t(\d+)g(\d+)v(\d+)modrt(\d+)b(\d+)')
    data_frame = []
    for i, fname in bar(enumerate(all_fnames)):
        #if i %100 == 0: print i,
        s = pattern.search(fname)
        mm, mh, cm, ch, om, alpha, teff, logg, micro, rot, res = s.groups()
        if mm == 'm':
            mscale = -1.0
        else:
            mscale = 1.0
        if cm == 'm':
            cscale = -1.0
        else:
            cscale = 1.0
        if om == 'm':
            oscale = -1.0
        else:
            oscale = 1.0

        teff = float(teff)
        logg = float(logg)
        mh = mscale * float(mh) / 10.
        ch = cscale * float(cscale) / 10.
        alpha = oscale * float(alpha) / 10.

        micro = float(micro)
        rot = float(rot)
        res = float(res)

        data_frame.append([teff, logg, mh, alpha, ch, rot, micro, res, fname])

    bosz_index.iloc[:, :] = data_frame

    return bosz_index


def make_grid_info(fname):
    """
    Make the HDF5 Grid Info file

    Parameters
    ----------
    fname: str

    """

    raw_index = make_raw_index()
    wavelength = np.loadtxt(raw_index.loc[0, 'filename'], usecols=(1,), unpack=True)

    with pd.HDFStore(fname) as fh:
        fh['index'] = raw_index
        fh['wavelength'] = pd.DataFrame(wavelength)
        fh['meta'] = pd.Series(bosz_meta)

def convert_bz2_memmap(fname):
    """
    Convert a bz2 file to memmap
    Parameters
    ----------
    fname : str

    Returns
    -------

    """
    fname_npy = fname.replace('.bz2', '.v1.npy')
    if os.path.exists(fname_npy):
        pass
    else:
        flux = pd.read_csv(fname, header=-1, usecols=(1,), delim_whitespace=True, dtype=np.float64)
        flux = flux.values[:,0]
        np.save(fname_npy, flux)

def cache_bosz_grid(delete=False):
    """
    Extract and cache BOSZ grid
    Parameters
    ----------
    delete: bool
        will delete the existing file

    Returns
    -------

    """
    all_fnames = glob('ascii/insbroad_300000/*/*/*/*')
    bar = ProgressBar(max_value=len(all_fnames))
    for i, fname in bar(enumerate(all_fnames)):
        convert_bz2_memmap(fname)

