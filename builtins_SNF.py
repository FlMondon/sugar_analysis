# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:34 2017

@author: mondon
"""

import sncosmo 
import astropy.units as u
import numpy as np
import pylab as P

def register_SNf_bands():
    
    wl = np.arange(3000,9500,0.001)
    transmission_U = np.zeros(len(wl))
    transmission_B = np.zeros(len(wl))
    transmission_V = np.zeros(len(wl))
    transmission_R = np.zeros(len(wl))
    transmission_I = np.zeros(len(wl))

    for i in range (len(wl)): 
        if wl[i] >= 3300 and wl[i] <= 4102:
            transmission_U[i]= 1.
        if wl[i] >= 4102 and wl[i] <= 5100:
            transmission_B[i]= 1. 
        if wl[i] >= 5200 and wl[i] <= 6289:
            transmission_V[i]= 1.
        if wl[i] >= 6289 and wl[i] <= 7607 :
            transmission_R[i]= 1.
        if wl[i] >= 7607 and wl[i] <= 9200 :
            transmission_I[i]= 1.        
            
    band_U = sncosmo.Bandpass(wl,transmission_U,wave_unit=u.AA,name='USNf')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl,transmission_B,wave_unit=u.AA,name='BSNf')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl,transmission_V,wave_unit=u.AA,name='VSNf')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl,transmission_R,wave_unit=u.AA,name='RSNf')    
    sncosmo.registry.register(band_R)
    band_I = sncosmo.Bandpass(wl,transmission_I,wave_unit=u.AA,name='ISNf')    
    sncosmo.registry.register(band_I)
    return wl, transmission_U, transmission_V, transmission_B, transmission_R, transmission_I
    
def plot_SNf_filters(wl, transmission_U, transmission_V, transmission_B, transmission_R, transmission_I):
    P.plot(wl, transmission_U)
    P.plot(wl, transmission_B)
    P.plot(wl, transmission_V)
    P.plot(wl, transmission_R)
    P.plot(wl, transmission_I)
    P.ylim(0,2)
    P.show()