import re
import os

from glob import glob

from astropy.io import fits
from astropy.table import Table
import pandas as pd
import numpy as np
from progressbar import ProgressBar

gotberg23_bibtex = """
@article{Gotberg:2023,
	abstract = {Massive stars (~8-25 M <SUB>⊙</SUB>) stripped of their hydrogen-rich envelopes via binary interaction are thought to be the main progenitors for merging neutron stars and stripped-envelope supernovae. We recently presented the discovery of the first set of such stripped stars in a companion paper. Here, we fit the spectra of 10 stars with new atmosphere models in order to constrain their stellar properties precisely. We find that the stellar properties align well with the theoretical expectations from binary evolution models for helium-core burning envelope-stripped stars. The fits confirm that the stars have high effective temperatures (T <SUB>eff</SUB> ~ 50-100 kK), high surface gravities ( $\mathrm log g\sim $ 5), and hydrogen-poor/helium-rich surfaces (X <SUB>H,surf</SUB> ~ 0-0.4) while showing for the first time a range of bolometric luminosities (10<SUP>3</SUP>-10<SUP>5</SUP> L <SUB>⊙</SUB>), small radii (~0.5-1 R <SUB>⊙</SUB>), and low Eddington factors (Γ<SUB> e </SUB> ~ 0.006-0.4). Using these properties, we derive intermediate current masses (~1-8 M <SUB>⊙</SUB>), which suggest that their progenitors were massive stars (~5-25 M <SUB>⊙</SUB>) and that a subset will reach core-collapse, leaving behind neutron stars or black holes. Using the model fits, we also estimate the emission rates of ionizing photons for these stars, which agree well with previous model expectations. Further, by computing models for a range of mass-loss rates, we find that the stellar winds are weaker than predicted by any existing scheme ( $ \dot M  _ \mathrm wind  \lesssim  10 ^ -9 $ M <SUB>⊙</SUB> yr<SUP>-1</SUP>). The properties of this first sample of intermediate-mass helium stars suggest they both contain progenitors of type Ib and IIb supernovae, and provide important benchmarks for binary evolution and population synthesis models.},
	adsnote = {Provided by the SAO/NASA Astrophysics Data System},
	adsurl = {https://ui.adsabs.harvard.edu/abs/2023ApJ...959..125G},
	archiveprefix = {arXiv},
	author = {{G{\"o}tberg}, Y. and {Drout}, M.~R. and {Ji}, A.~P. and {Groh}, J.~H. and {Ludwig}, B.~A. and {Crowther}, P.~A. and {Smith}, N. and {de Koter}, A. and {de Mink}, S.~E.},
	date-added = {2024-01-17 16:02:23 -0800},
	date-modified = {2024-01-17 16:02:23 -0800},
	doi = {10.3847/1538-4357/ace5a3},
	eid = {125},
	eprint = {2307.00074},
	journal = {\apj},
	keywords = {Binary stars, Close binary stars, Interacting binary stars, Early-type stars, Helium-rich stars, Helium burning, Stellar properties, Stellar spectral types, Stellar spectral lines, Ionization, Stellar winds, 154, 254, 801, 430, 715, 716, 1624, 2051, 1630, 2068, 1636, Astrophysics - Solar and Stellar Astrophysics, Astrophysics - Astrophysics of Galaxies},
	month = dec,
	number = {2},
	pages = {125},
	primaryclass = {astro-ph.SR},
	title = {{Stellar Properties of Observed Stars Stripped in Binaries in the Magellanic Clouds}},
	volume = {959},
	year = 2023,
	bdsk-file-1 = {YnBsaXN0MDDSAQIDBFxyZWxhdGl2ZVBhdGhZYWxpYXNEYXRhbxAaAEcAbwMIAHQAYgBlAHIAZwAvAEcAbwMIAHQAYgBlAHIAZwBfADIAMAAyADMALgBwAGQAZk8RAdIAAAAAAdIAAgAADE1hY2ludG9zaCBIRAAAAAAAAAAAAAAAAAAAAOHGKAxCRAAB/////xBHmnRiZXJnXzIwMjMucGRmAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAD/////4c2qjwAAAAAAAAAAAAEAAwAACiBjdQAAAAAAAAAAAAAAAAAHR5p0YmVyZwAAAgBpLzpVc2VyczphYmhpbWF0OkxpYnJhcnk6TW9iaWxlIERvY3VtZW50czpjb21+YXBwbGV+Q2xvdWREb2NzOkFjYWRlbWljOlBhcGVyczpHb8yIdGJlcmc6R2/MiHRiZXJnXzIwMjMucGRmAAAOACQAEQBHAG8DCAB0AGIAZQByAGcAXwAyADAAMgAzAC4AcABkAGYADwAaAAwATQBhAGMAaQBuAHQAbwBzAGgAIABIAEQAEgBnVXNlcnMvYWJoaW1hdC9MaWJyYXJ5L01vYmlsZSBEb2N1bWVudHMvY29tfmFwcGxlfkNsb3VkRG9jcy9BY2FkZW1pYy9QYXBlcnMvR2/MiHRiZXJnL0dvzIh0YmVyZ18yMDIzLnBkZgAAEwABLwAAFQACAA7//wAAAAgADQAaACQAWwAAAAAAAAIBAAAAAAAAAAUAAAAAAAAAAAAAAAAAAAIx},
	bdsk-url-1 = {https://doi.org/10.3847/1538-4357/ace5a3},
	bdsk-url-2 = {https://ui.adsabs.harvard.edu/abs/2023ApJ...959..125G},
	bdsk-url-3 = {https://ui.adsabs.harvard.edu/link_gateway/2023ApJ...959..125G/EPRINT_HTML}
}

@dataset{Gotberg:2023a,
	author = {G{\"o}tberg, Y. and Drout, M.R. and Ji, A.P. and Groh, J.H. and Ludwig, B.A. and Crowther, P.A. and Smith, N. and de Koter, A. and de Mink, S.E.},
	date-added = {2024-01-22 09:12:21 -0800},
	date-modified = {2024-01-22 09:12:39 -0800},
	doi = {10.5281/zenodo.7976200},
	month = jul,
	publisher = {Zenodo},
	title = {{Stellar properties of observed stars stripped in binaries in the Magellanic Clouds - Spectral Models}},
	url = {https://doi.org/10.5281/zenodo.7976200},
	version = 1,
	year = 2023,
	bdsk-url-1 = {https://doi.org/10.5281/zenodo.7976200}
}

"""

gotberg23_meta = {
    'bibtex':gotberg23_bibtex,
    'parameters':[
        'teff', 'logg',                     # Teff and logg at tau=2/3
        'teff_tau20', 'logg_tau20',         # Teff and logg at tau=20
        'xh_surf', 'xhe_surf', 'z_surf',    # Hydrogen, Helium, and Z abundance
    ],
    'wavelength_unit':'Angstrom',
    'wavelength_type':'vacuum',
    'flux_unit': 'erg/s/cm^2/angstrom',
}


def make_raw_index(
        res=300000.0,
        grid_params_file='S41_spectral_model_grid_parameters.txt',
    ):
    """
    Read all Gotberg+ 23 grid files and generate a raw index with filename association.

    Returns
    -------
        bosz_index : pd.DataFrame
    """
    
    # Read in spectral model grid parameters file as Table
    params_table = Table.read(
        grid_params_file,
        format='ascii',
        delimiter='\s',
        data_start=2,
    )
    
    
    # Construct file names from model name files
    models_directory = 'S41_spectra_spectral_model_grid'
    
    filenames = list(params_table['Model_name'])
    
    for i, model_name in enumerate(params_table['Model_name']):
        filenames[i] = '{0}/SED_{1}.txt'.format(
            models_directory,
            model_name.replace('Teff', 'T'),
        )
    
    # Extract all other parameters from the params table
    teff_arr = params_table['Teff(tau=2/3)'].data
    logg_arr = params_table['log_geff(tau=2/3)'].data
    rad_arr = params_table['Reff(tau=2/3)'].data
    
    teff_tau20_arr = params_table['Teff(tau=20)'].data
    logg_tau20_arr = params_table['log_geff(tau=20)'].data
    rad_tau20_arr = params_table['Reff(tau=20)'].data
    
    xh_surf_arr = params_table['X_H'].data
    xhe_surf_arr = params_table['X_He'].data
    z_surf_arr = params_table['Z'].data
    
    res_arr = np.full_like(teff_arr, fill_value=res)
    
    # Construct and output pandas data frame
    output_pd_df = pd.DataFrame({
        'teff':teff_arr,
        'logg':logg_arr,
        'rad':rad_arr,
        'teff_tau20':teff_tau20_arr,
        'logg_tau20':logg_tau20_arr,
        'rad_tau20':rad_tau20_arr,
        'xh_surf':xh_surf_arr,
        'xhe_surf':xhe_surf_arr,
        'z_surf':z_surf_arr,
        'res':res_arr,
        'filename':filenames,
    })
    
    return output_pd_df

def make_grid_info(fname='./gotberg23_info_v1p0.h5'):
    """
    Make the HDF5 Grid Info file
    
    Parameters
    ----------
    fname: str
    
    """
    
    raw_index = make_raw_index()
    
    # Use the first model file's wavelength column as reference
    
    wtab = Table.read(
        raw_index.loc[0, 'filename'],
        format='ascii.commented_header',
    )
    
    wavelength = wtab['Wavelength']
    
    with pd.HDFStore(fname) as fh:
        fh['index'] = raw_index
        fh['wavelength'] = pd.DataFrame(wavelength)
        fh['meta'] = pd.Series(gotberg23_meta)

def convert_sed_memmap(fname):
    """
    Convert a Gotberg23 SED file to memmap
    Parameters
    ----------
    fname : str

    Returns
    -------

    """
    fname_npy = fname.replace('.txt', '.v1.npy')
    if os.path.exists(fname_npy):
        pass
    else:
        sed_tab = Table.read(
            fname,
            format='ascii.commented_header',
        )
        
        flux = sed_tab['Flambda'].value
        
        np.save(fname_npy, flux)

def cache_gotberg23_grid(delete=False):
    """
    Extract and cache Gotberg23 grid
    Parameters
    ----------
    delete: bool
        will delete the existing file

    Returns
    -------

    """
    all_fnames = glob('S41_spectra_spectral_model_grid/SED*.txt')
    bar = ProgressBar(maxval=len(all_fnames))
    for i, fname in bar(enumerate(all_fnames)):
        convert_sed_memmap(fname)
