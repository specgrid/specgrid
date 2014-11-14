# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
This is an Astropy affiliated package.
"""

# Affiliated packages may add whatever they like to this file, but
# should keep this content at the top.
# ----------------------------------------------------------------------------
from ._astropy_init import *
# ----------------------------------------------------------------------------

# For egg_info test builds to pass, put package imports here.
if not _ASTROPY_SETUP_:
    from specgrid.fix_spectrum1d import Spectrum1D
    from specgrid.base import SpectralGrid
    from specgrid import plugins
    from specgrid.model_star import (ModelStar, ModelInstrument, Observation,
                                     assemble_observation)


import logging

logger = logging.getLogger('specgrid')
logger.setLevel(logging.INFO)
console_handler = logging.StreamHandler()
console_formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(console_formatter)
logger.addHandler(console_handler)