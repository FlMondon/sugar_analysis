#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 10:40:00 2019

@author: florian
"""

import sncosmo
import numpy as np
import os
from iminuit import Minuit
import pickle as pkl
from .constant import t_min_sug, t_max_sug 
from .fitting_lc import LC_Fitter
from .Hubble_fit import read_input_data_SNf 
from .load_sugext import register_SUGAR_ext

class train_transext(object):
    
    def __init__(self, modelcov=True, sad_path = '../../',
                 modeldir='../../sugar_model/', mod_errfile='model_err_sug.dat',
                 width=10, param_sug =['q1', 'q3']): 
    
            register_SUGAR_ext(modeldir=modeldir, mod_errfile=mod_errfile,
                       version='0.0')

    