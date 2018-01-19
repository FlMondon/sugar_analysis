# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:34 2017

@author: mondon
"""

import sncosmo 
import astropy.units as u
import numpy as np
import pylab as P
import pyfits
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d

def register_SNf_bands():
    
    wl = np.arange(3000,9500,0.1)
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

def load_spectral_magsys_fits2(relpath, name=None, version=None):
    fits = pyfits.open(relpath)
    fit = fits[1]
    wl = np.zeros(len(fit.data))
    flux = np.zeros(len(fit.data))
    for i in range(len(fit.data)):
         wl[i] = fit.data[i][0]
         flux[i] = fit.data[i][1]
    refspectrum = sncosmo.spectrum.Spectrum(wl, flux,
                           unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    return sncosmo.magsystems.SpectralMagSystem(refspectrum, name=name)
    

def mag_sys_SNF():
    """
    define magnitude systeme for snf
    """
    
    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0', load_spectral_magsys_fits2,
                             args=['../sugar_analysis_data/Vega.fits'],
                             meta={'description': 'use vega spectrum that come from snf data'})        
    bands_snf ={'USNF': ('vega_snf_0', 0.),
	'BSNF': ('vega_snf_0', 0.),
	'VSNF': ('vega_snf_0', 0.),
	'RSNF': ('vega_snf_0', 0.),
	'ISNF': ('vega_snf_0', 0)}
    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_snf),'vega_snf') 
 
def plot_SNf_filters(wl, transmission_U, transmission_V, transmission_B, transmission_R, transmission_I):
    P.plot(wl, transmission_U)
    P.plot(wl, transmission_B)
    P.plot(wl, transmission_V)
    P.plot(wl, transmission_R)
    P.plot(wl, transmission_I)
    P.ylim(0,2)
    P.show()

   
def define_zp_snf(wl,transmission, path_ref_spectrum):
    """
    """
    wl_min = 3000
    wl_max = 9500
    fits = pyfits.open(path_ref_spectrum)
    fit = fits[1]
    wlref = np.zeros(len(fit.data))
    fluxref = np.zeros(len(fit.data))
    for i in range(len(fit.data)):
         wlref[i] = fit.data[i][0]
         fluxref[i] = fit.data[i][1]
    
    splref = Spline1d(wlref, fluxref, k=1,ext = 1)
    splband = Spline1d(wl, transmission, k=1,ext = 1)

    #computation of the integral
    dt = 100000
    xs = np.linspace(float(wl_min), float(wl_max), dt)
    dxs = (float(wl_max-wl_min)/(dt-1))    
    
    I = np.sum(splref(xs)*xs*splband(xs)*dxs)
    
    zp = 2.5*np.log10(1/I)
    return zp 