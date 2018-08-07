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
from .builtins_SNF import register_SNf_bands_width, mag_sys_SNF_width, register_SUGAR
from .read_data import read_meta_SNF
