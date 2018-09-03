#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 19:34:46 2018

@author: florian
"""

""" test the photometric fit. """

import numpy as np
from sugar_analysis import sugar_generator as sg
sugar_model = '/home/florian/sugar_model/'

class test_fit(object):
    """
    unit test for photometric light curve fit and spectroscopic light curve fit
    """
    def __init__(self, q1, q2, q3, Av, grey):
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3
        self.Av = Av
        self.grey = grey
        self.ss = sg.sugar_spectrum()
        self.ss._parameters = np.array([q1, q2, q3, Av, grey])
        self.model =  np.loadtxt(sugar_model+'SUGAR_model_v1.asci')

        
    def test_fit_phot(self):
        self.data_phot = self.ss.AstropyTable_flux()
        res_phot, fitted_model = self.ss.fit_lc_sugar(self.data_phot)
        res = np.array([res_phot.parameters[3], 
                  res_phot.parameters[4], 
                  res_phot.parameters[5], 
                  res_phot.parameters[6], 
                  -2.5*np.log10(res_phot.parameters[2])])
        np.testing.assert_allclose(res, self.ss._parameters, atol=1e-3)

    def test_fit_spec(self):
        h = self.ss.fit_spec_sugar()
        np.testing.assert_allclose(h, self.ss._parameters, atol=1e-3)
        
if __name__=="__main__":
    tf = test_fit(1.,1.,1.,0.1,38.)
    tf.test_fit_phot()
    tf.test_fit_spec()
    