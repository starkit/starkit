import re
import os

from glob import glob

from astropy.io import fits
import pandas as pd
import numpy as np

from starkit.gridkit.io.pollux.base import make_grid_info,cache_pollux_grid


def main():
    spectra_dir = '/Volumes/data/pollux_grids/cmfgen_martins/IR'

    make_grid_info('pollux_cmfgen_info_IR.h5',spectra_dir=spectra_dir)
    cache_pollux_grid(spectra_dir=spectra_dir)


if __name__ == "__main__":
    main()