# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an Astropy affiliated package.
"""

import os, sys
on_rtd = os.environ.get('READTHEDOCS') == 'True'
from ._astropy_init import *

if not on_rtd:
    # Affiliated packages may add whatever they like to this file, but
    # should keep this content at the top.
    # ----------------------------------------------------------------------------
    
    # ----------------------------------------------------------------------------

    # For egg_info test builds to pass, put package imports here.
    from starkit.utils.colored_logger import ColoredFormatter, formatter_message
    from starkit.base.assemble_model import assemble_model
    from starkit.fix_spectrum1d import SKSpectrum1D
    import starkit.base.operations as operations

    import logging

    FORMAT = "[$BOLD%(name)-20s$RESET][%(levelname)-18s]  %(message)s ($BOLD%(filename)s$RESET:%(lineno)d)"
    COLOR_FORMAT = formatter_message(FORMAT, True)

    logging.captureWarnings(True)
    logger = logging.getLogger('starkit')
    logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler(sys.stdout)
    console_formatter = ColoredFormatter(COLOR_FORMAT)
    console_handler.setFormatter(console_formatter)

    logger.addHandler(console_handler)
    logging.getLogger('py.warnings').addHandler(console_handler)
