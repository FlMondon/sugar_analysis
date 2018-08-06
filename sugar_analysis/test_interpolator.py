#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:54:41 2018

@author: florian
"""
import numpy as np
import sncosmo
import cPickle

CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27        # [erg s]
p2 = '../sugar_model/'

#                    
            
def _flux(phase, wave):


    mod_flux = cPickle.load(open(p2+'mod_flux.pkl'))
    flux_inter = sncosmo.salt2utils.BicubicInterpolator(phase, wave, mod_flux)

    
    return flux_inter(phase,wave)