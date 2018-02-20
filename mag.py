#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 21:56:02 2018

@author: florian
"""

import numpy as np
from math import log, log10
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
import pyfits

CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27       # [erg s]
# wavelength limits for salt2 model
wl_min_sal = 3000
wl_max_sal = 7000
spec_path = '../sugar_analysis_data/SNF-0203-CABALLO2/validation/CSS110918_01/CSS110918_01_P002391_merged.fits'

class Spectrum():
    """
    """
    def __init__(self, photons=True):
        """
        """
        self.photons = photons
        
    def read_spectrum(self):
        """
        """
        
        self._fits = pyfits.open(spec_path)
        self.spec = self._fits[0]
        self._hdr = self.spec.header.copy()       # Spectrum header

        self._hdr['CRPIX1'] = self._hdr.get('CRPIX1', 1)  # Make it mandatory

        self.npts = self._hdr['NAXIS1']
        self.step = self._hdr['CDELT1']
        self.start= self._hdr['CRVAL1'] - (self._hdr['CRPIX1'] - 1) * self.step
        self.end  = self.start + (self.npts - 1) * self.step
        self.x = np.linspace(self.start, self.end, self.npts)  # Wavelength
        self.y = self.spec.data.copy()                         # Signal

        self.splspec = Spline1d(self.x, self.y, k=1,ext = 1)    
    
    
    def flux(self,band):
        """
        """
    
        self.read_spectrum()
        if band == 'USNf' :
            band_file = 'USNf_3300-4102.dat'
            
        if band == 'BSNf' :
            band_file = 'BSNf_4102-5100.dat'
            
        if band == 'VSNf' :
            band_file = 'VSNf_5200-6289.dat'
            
        if band == 'RSNf' :
            band_file = 'RSNf_6289-7607.dat'
            
        if band == 'ISNf' :
            band_file = 'ISNf_7607-9200.dat'
            
            
        filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file)
        wlen = filt2[:,0]
        tran = filt2[:,1]
        self.splB = Spline1d(wlen, tran, k=1,ext = 1)    
        
        #computation of the integral
        dt = 100000
        xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
        dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))
    
        self.flux = np.sum(self.splspec(xs)*xs/ (HPLANCK * CLIGHT)*self.splB(xs)*dxs)
#        print self.flux
    
    
    def zero_point (self, band):
        """
        """
        self.read_spectrum()
        if band == 'USNf' :
            lmax = 4102.0
            lmin = 3300.0
            
        if band == 'BSNf' :
            lmax = 5100.0
            lmin = 4102.0 
            
        if band == 'VSNf' :
            lmax = 6289.0
            lmin = 5200.0
            
        if band == 'RSNf' :
            lmax = 7607.0
            lmin = 6289.0
            
        if band == 'ISNf' :
             lmax = 9200.0
             lmin = 7607.0
            

    
        self.zp = log(lmax / lmin) / HPLANCK
        self.zp = 2.5*log10(self.zp)
#        print self.zp
        


        
    def _mref(self,band):
        fits = pyfits.open('../sugar_analysis_data/Vega.fits')
        fit = fits[1]
        wl = np.zeros(len(fit.data))
        flux = np.zeros(len(fit.data))
        for i in range(len(fit.data)):
             wl[i] = fit.data[i][0]
             flux[i] = fit.data[i][1]
            
        
        self.splref = Spline1d(wl, flux, k=1,ext = 1)
        
        #computation of the integral
        dt = 100000
        xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
        dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))
        
        self.zero_point(band)
        self.I2 = np.sum(self.splref(xs)*xs/ (HPLANCK * CLIGHT)*self.splB(xs)*dxs)       
        print self.I2
        print 2.5*log10(self.I2)
        self.mref = -2.5*log10(self.I2) + self.zp - 48.59 
#        print self.mref
        
        
        
    def flux_to_mag(self,band):
        """
        """
        self.flux(band)
#        if band == 'USNf' :
#            mref = 9.787
#        if band == 'BSNf' :
#            mref = 9.791
#        if band == 'VSNf' :
#            mref = 9.353
#        if band == 'RSNf' :
#            mref = 9.011
#        if band == 'ISNf' :
#            mref = 8.768
        #zeropint
        self.zero_point(band)
        self._mref(band)


        self.m = -2.5*log10(self.flux) + self.zp - 48.59 - self.mref 
#        print -2.5*log10(self.flux)
        print self.m
        print self.flux
        
    
    
    
        
        
        
        
        