import re
import os

from glob import glob

from astropy.io import fits
import pandas as pd
import numpy as np
from progressbar import ProgressBar
from .parse_filename import parse_model_filename

import json
import io # Use the io module for cross-compatible file handling

pollux_bibtex = """
@ARTICLE{pollux REFERENCE HERE
}

"""

'''
Meta data includes
{
    "metallicity": 0.0, 
    "spectral_domain": "VIS", 
    "t_eff": 35000, 
    "n_fe": 0.0, 
    "gravity": 3.4, 
    "model_type": "s", 
    "c_fe": 0.0, 
    "o_fe": 0.0, 
    "beta": 0.8, 
    "terminal_velocity": 1760, 
    "alpha_fe": 0.0, 
    "model_family": "C", 
    "mass_loss": -5.41, 
    "microturbulent_velocity": 5.0, 
    "first_clumping_parameter": 1, 
    "second_clumping_parameter": 0
}
'''

meta_example = {
    "metallicity": 0.0, 
    "spectral_domain": "VIS", 
    "t_eff": 35000, 
    "n_fe": 0.0, 
    "gravity": 3.4, 
    "model_type": "s", 
    "c_fe": 0.0, 
    "o_fe": 0.0, 
    "beta": 0.8, 
    "terminal_velocity": 1760, 
    "alpha_fe": 0.0, 
    "model_family": "C", 
    "mass_loss": -5.41, 
    "microturbulent_velocity": 5.0, 
    "first_clumping_parameter": 1, 
    "second_clumping_parameter": 0
}

pollux_meta = {'bibtex':pollux_bibtex,
             'parameters':['teff','logg','mh'], # we will limit the number of parameters to make the interpolator faster
             'wavelength_unit':'Angstrom',
             'wavelength_type':'vacuum',
             'flux_unit': 'erg/s/cm^2/angstrom'}


def make_raw_index(mh=-0.08,alpha=-0.05,res=300000.0,spectra_dir='/Volumes/data/pollux_grids/cmfgen_martins/IR'):
    """
    Read all Pollux files and generate a raw index with filename association. This is targetting the CMFGEN Pollux grid.

    Returns
    -------
        pollux_index : pd.DataFrame
    """
    all_fnames = glob(os.path.join(spectra_dir,'*.spec'))

    nfiles = len(all_fnames)
    mh_arr = np.zeros(nfiles)
    alpha_arr = np.zeros(nfiles)
    teff_arr = np.zeros(nfiles)
    logg_arr = np.zeros(nfiles)
    micro_arr = np.zeros(nfiles)
    nfe_arr = np.zeros(nfiles)
    cfe_arr = np.zeros(nfiles)
    ofe_arr = np.zeros(nfiles)
    vturb_arr = np.zeros(nfiles)
    vinfty_arr = np.zeros(nfiles)
    mass_arr = np.zeros(nfiles)
    lum_arr = np.zeros(nfiles)

    for i in np.arange(nfiles):
        filename = all_fnames[i]
        

        p = filename + '.txt'
        params = parse_cmfgen_meta_file(p)
        
        if params is None:
            continue
        mh = params['metallic']
        alpha = params['alpha']
        mdot = params['Mdot']
        c = params['carbon']
        n = params['nitrogen']
        o = params['oxygen'] 
        micro = params['turbvel_Vmin'] # microturbulent velocity (min) km/s
        teff = params['Teff'] # effective temperature (K)
        logg = params['logg'] # log10 cm/s^2
        mass = params['mass'] # mass (solar masses)
        lum = params['lum'] # luminosity (solar luminosities)
        vinfty = params['Vinfty'] # terminal velocity

        #he_h = params['he_h']
        mh_arr[i] = mh
        alpha_arr[i] = alpha
        teff_arr[i] = teff
        logg_arr[i] = logg
        micro_arr[i] = micro
        nfe_arr[i] = n
        cfe_arr[i] = c
        ofe_arr[i] = o
        vinfty_arr[i] = vinfty
        mass_arr[i] = mass
        lum_arr[i] = lum

    return pd.DataFrame({'filename':all_fnames,
                         'mh':mh_arr,
                         'teff':teff_arr,
                         'logg':logg_arr})




def parse_cmfgen_meta_file(filepath):
    """
    Opens, reads, and parses a CMFGEN spec file into a Python dictionary.
    This function is compatible with both Python 2 and 3.

    Args:
        filepath (str): The path to the .spec.txt file.

    Returns:
        dict: A dictionary with the parsed key-value pairs, or an empty
        dictionary if the file cannot be found.
    """
    try:
        # Use io.open for consistent encoding support across Python 2 and 3
        with io.open(filepath, 'r', encoding='utf-8') as f:
            file_content = f.read()
    except IOError:
        # Catch IOError, which handles file-not-found errors in both versions.
        # Use .format() for string formatting to maintain compatibility.
        print("Error: The file at path '{0}' was not found.".format(filepath))
        return {} # Return an empty dictionary on error

    parsed_data = {}
    # Regex to capture key = 'value' pairs
    pattern = re.compile(r"^\s*([\w_]+)\s*=\s*'([^']*)'")

    # Iterate over each line in the file content
    for line in file_content.splitlines():
        match = pattern.match(line)
        if match:
            key = match.group(1)
            value = match.group(2)
            parsed_data[key] = value
            
    return parsed_data

def make_grid_info(fname,spectra_dir='/Volumes/data/pollux_grids/cmfgen_martins/IR'):
    """
    Make the HDF5 Grid Info file

    Parameters
    ----------
    fname: str

    """

    raw_index = make_raw_index(spectra_dir=spectra_dir)
    wtab = pd.read_csv(
        raw_index.loc[0, 'filename'],
        sep=r'\s+',
        header=None,
        names=['wavelength', 'flux'],
        usecols=(0, 1),
    )
    wavelength = wtab['wavelength']

    with pd.HDFStore(fname) as fh:
        fh['index'] = raw_index
        fh['wavelength'] = pd.DataFrame(wavelength)
        fh['meta'] = pd.Series(pollux_meta)

def convert_bz2_memmap(fname):
    """
    Convert a bz2 file to memmap
    Parameters
    ----------
    fname : str

    Returns
    -------

    """
    fname_npy = fname.replace('.spec', '.spec.v1.npy')
    if os.path.exists(fname_npy):
        pass
    else:
        flux = pd.read_csv(
            fname,
            usecols=(1,),
            sep=r'\s+',
            dtype=np.float64,
            header=None,
            names=['flux'],
        )
        flux = flux.values[:,0]
        np.save(fname_npy, flux)

def cache_pollux_grid(delete=False,spectra_dir='/Volumes/data/pollux_grids/cmfgen_martins/IR'):
    """
    Extract and cache pollux grid
    Parameters
    ----------
    delete: bool
        will delete the existing file

    Returns
    -------

    """
    all_fnames = glob(os.path.join(spectra_dir,'*.spec'))
    bar = ProgressBar(maxval=len(all_fnames))
    for i, fname in bar(enumerate(all_fnames)):
        convert_bz2_memmap(fname)
