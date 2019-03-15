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
from .cosmo_tools import distance_modulus_th, int_cosmo, luminosity_distance
from .builtins import register_SNf_bands_width, mag_sys_SNF_width,  builtins_jla_bandpasses, mag_sys_jla
from .load_sugar import register_SUGAR
from .Hubble_fit import read_input_data_SNf, get_hubblefit, generate_fake_data
from .fitting_lc import LC_Fitter
from .sugar_generator import sugar_simulation
from .data_table import build_data
from .plot_lc import plot_lc_res
