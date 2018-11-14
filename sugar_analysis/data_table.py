#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:55:09 2018

@author: florian
"""


import numpy as np
import sncosmo
from sugar_analysis import constant as cst
from astropy.table import Table
#from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from sugar_analysis import builtins
import cPickle as pkl
output_path = '../../'

    
class build_data(object):

    def __init__(self,meta=output_path+'sugar_analysis_data/SNF-0203-CABALLO2/META.pkl'):
        builtins.register_SNf_bands_width(width=10)
        builtins.mag_sys_SNF_width(width=10)
        builtins.register_SUGAR()  
        self.noTH = True
        self.spec = pkl.load(open(output_path+
                                  'sugar_analysis_data/SNF-0203-CABALLO2/meta_spectra.pkl'))
        self.meta= pkl.load(open(meta))
      
    def integral_to_phot(self, wave, flux, var=None):

        self.start = min(wave)
        self.end =  max(wave)
        # Interpolate filter over spectrum (natively regularly sampled)
        f = self.interpolate(wave)
        f *= wave/(cst.CLIGHT*1.0e13*cst.HPLANCK)
        
        #computation of the integral
        dxs = (float(self.end - self.start)/(len(wave)-1))
        phot = np.sum(flux*f*dxs)
        if  isinstance(var, np.ndarray):
            dphot = np.sqrt((var*f**2).sum())*dxs
        else:
            dphot = None
#        for i in range(len(f)):
#            if f[i] > 1:
#                print '%e'%(f[i])
#        print '############################################'
        return phot, dphot
    
    def interpolate(self, x):
        """
        Interpolate filter over input array x (not necessarily
        regularly sampled).

        Beware that interpolation doesn't go wild *outside* the filter
        definition range, or even *inside*.
        """

        from scipy.interpolate import UnivariateSpline
        filt2 = np.genfromtxt(output_path+
                              'sugar_analysis_data/data/Instruments/Florian/'+self.band+'.dat')
        transp = 0.000
        for i, trans in enumerate(filt2[:,1]):
            if trans > 0.000 and transp == 0.000:
                start = filt2[:,0][i]
            if transp > 0.000 and trans == 0.000:
                end = filt2[:,0][i-1]
            transp = trans
        y = self.pixel_fracs(x, start, end)
        
        if self.noTH :
            inside = y > 0
            spline = UnivariateSpline(filt2[:,0], filt2[:,1], s=0)
            y[inside] = np.maximum(spline(x[inside]), 0)
        return y

    def pixel_fracs(self, x, xmin, xmax):
        """Return which fraction of each (Voronoi) pixel centered on (not
        necessaryly regularly sampled) x is contained within wavelength
        range [xmin,xmax]."""
    
        x = np.asarray(x)        # Center of pixels
        n = len(x)
    
        # Pixel boundaries
        dx = np.diff(x) / 2        # Half px width
        xb = np.concatenate(([x[0] - dx[0]],    # 1st px left boundary
                            x[1:] - dx,        # Middle pixel left boundaries
                            [x[-1] + dx[-1]]))  # Last px right boundary
    
        # Fractions
        f = np.zeros(n)
        f[(xmin <= x) & (x <= xmax)] = 1.  # Rough approach 1st
    
        # Compute precise fractions at domain boundaries
        imin = xb.searchsorted(xmin)  # xb[imin-1] < xmin < xb[imin]
        imax = xb.searchsorted(xmax)  # xb[imax-1] < xmax < xb[imax]
        if 0 < imin <= n:
            f[imin - 1] = (min(xb[imin], xmax) - xmin) / (xb[imin] - xb[imin - 1])
        if 0 < imax <= n:
            f[imax - 1] = (xmax - max(xb[imax - 1], xmin)) / \
                (xb[imax] - xb[imax - 1])
#        if self.band == 'RSNf':
#            for i in range(len(f)):
#                if f[i]>0.1:
#                    print x[i], f[i]       
        return f

    def build_Astropy_Table(self, sn_name, band_used=['new_fU_10','fB_10','fV_10','fR_10','new_fI_10']):
        """
        """    
        time = []
        flux = []
        flux_err = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('vega_snf_10')

        
        for b in band_used:
            self.band = b
            if b == 'fB_10' or b == 'new_fU_10':
                spectra = 'spectra_B'
            elif b == 'fV_10' or  b == 'fR_10' or  b == 'new_fI_10':
                spectra = 'spectra_R'
            for p in self.spec[sn_name][spectra].keys():
                phot_flux, dphot = self.integral_to_phot(
                        self.spec[sn_name][spectra][str(p)]['wave'],
                        self.spec[sn_name][spectra][str(p)]['flux'],
                        var=self.spec[sn_name][spectra][str(p)]['var']) 
                flux.append(phot_flux)
                time.append(float(p))
                band.append(self.band)
                flux_err.append(dphot*self.meta[sn_name]['target.errorscale'])
                zp.append(2.5*np.log10(vega.zpbandflux(self.band)))
                zpsys.append('vega_snf_10')
                
        data = Table([time, band, flux, flux_err, zp, zpsys], 
                     names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), 
                     meta={'name': 'data'})
        return data