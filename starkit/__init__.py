# Licensed under a 3-clause BSD style license - see licenses/LICENSE.rst

"""
StarKit: a collection of tools for working with stellar data.
"""

import logging
import sys

from starkit.version import version as __version__, githash as __githash__
from starkit.utils.colored_logger import ColoredFormatter, formatter_message
from starkit.base.assemble_model import assemble_model
from starkit.fix_spectrum1d import SKSpectrum1D
import starkit.base.operations as operations

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
