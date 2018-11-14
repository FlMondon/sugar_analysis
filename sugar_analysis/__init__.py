#!/usr/bin/env python

"""
Some description.
"""

#import os
#import glob

# Automatically import all modules (python files)
#__all__ = [os.path.basename(m).replace('.py', '') for m in glob.glob("sugar/*.py")
#           if '__init__' not in m]

from .math_toolbox import comp_rms
from .builtins import register_SNf_bands_width, mag_sys_SNF_width, register_SUGAR,  builtins_jla_bandpasses, mag_sys_jla
from .Hubble_fit import read_input_data_SNf, get_hubblefit
from .fitting_lc import LC_Fitter
from .cosmo_tools import distance_modulus_th, int_cosmo, luminosity_distance
from .sugar_generator import sugar_simulation
from .data_table import build_data
