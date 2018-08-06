#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 31 23:35:31 2018

@author: florian
"""

import numpy as np 
#import pylab as p 
import lmfit as fit

class chisq():
    
    def __init__(self, data_x, data_y, data_err):
        self.data_y = data_y
        self.data_x = data_x
        self.data_err = data_err
        print 'toto'
        
    def model(self,params):
        return params[0]*self.data_x + params[1]
    
    def chisq(self,params):
        
        residuals = (self.data_y-self.model(params))/self.data_err
        return np.sum(residuals**2)
    
    def chisq_mcmc(self, params):
        return -self.chisq(params)
#    def fitting(self):
#        params = fit.Parameters()
#        params.add('a', value=0.)
#        params.add('b', value=0.)
#                
#        return fit.minimize(self.chisq,params)