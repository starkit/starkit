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
             'wavelength_type':'vacuum',
             'flux_unit': 'erg/s/cm^2/angstrom'}


def make_raw_index():
    """
    Read all Phoenix files and generate a raw index with filename association.

    Returns
    -------
        bosz_index : pd.DataFrame
    """
    all_fnames = glob('ascii/insbroad_300000/*/*/*/*.asc.bz2')

    nfiles = len(all_fnames)
    mh_arr = np.zeros(nfiles)
    ch_arr = np.zeros(nfiles)
    alpha_arr = np.zeros(nfiles)
    teff_arr = np.zeros(nfiles)
    logg_arr = np.zeros(nfiles)
    micro_arr = np.zeros(nfiles)
    rot_arr = np.zeros(nfiles)
    res_arr = np.zeros(nfiles)
    pattern = re.compile('a(mp|mm)(\d+)(cp|cm)(\d+)+(op|om)(\d+)t(\d+)g(\d+)v(\d+)modrt(\d+)b(\d+)')
    pattern_dir = re.compile('metal_(.....)\/carbon_(.....)\/alpha_(.....)')

    for i in np.arange(nfiles):
        filename = all_fnames[i]
        base = filename.split('.')[0]
        s = pattern.search(filename)
        s2 = pattern_dir.search(filename)
        mm,mh,cm,ch,om,alpha,teff,logg,micro,rot,res =  s.group(1,2,3,4,5,6,7,8,9,10,11)
        mh,ch,alpha = s2.group(1,2,3) # use the directory names for more accurate grid points

        logg = logg[0]+'.'+logg[1:]
        micro = micro[0]+'.'+micro[1:]

        mh_arr[i] = float(mh)
        ch_arr[i] = float(ch)
        alpha_arr[i] = float(alpha)
        teff_arr[i] = float(teff)
        logg_arr[i] = float(logg)
        micro_arr[i] = float(micro)
        rot_arr[i] = float(rot)
        res_arr[i] = float(res)

    return pd.DataFrame({'mh':mh_arr,'ch':ch_arr,'alpha':alpha_arr,'teff':teff_arr,
                      'logg':logg_arr,'micro':micro_arr,'rot':rot_arr,
                      'res':res_arr,'filename':all_fnames})

def make_grid_info(fname):
    """
    Make the HDF5 Grid Info file

    Parameters
    ----------
    fname: str

    """

    raw_index = make_raw_index()
    wavelength = np.loadtxt(raw_index.loc[0, 'filename'], usecols=(0,), unpack=True)

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
