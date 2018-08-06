##!/usr/bin/env python2
## -*- coding: utf-8 -*-
#"""
#Created on Wed May 23 14:36:30 2018
#
#@author: florian
#"""
#
#import numpy as np 
#import sncosmo 
#
##model path
#p2 = '../sugar_model/'
#
#CLIGHT = 2.99792458e18         # [A/s]
#HPLANCK = 6.62606896e-27       # [erg s]
## wavelength limits for sugar model
#wl_min_sug = 3254.01639
#wl_max_sug = 8649.03871
#wl_min_salt = 2000.000000
#wl_max_salt = 9200.000000
#
#class spectrum_simulator():
#
#    def __init__(self):
#        """
#        """
#        self.m0file='sugar_template_0.dat'
#        self.m1file='sugar_template_1.dat'
#        self.m2file='sugar_template_2.dat'
#        self.m3file='sugar_template_3.dat'
#        self.m4file='sugar_template_4.dat'        
#        self._parameters = np.array([0., 0., 0., 0., 37.])
#        self._SCALE_FACTOR = 1.
#        self._model = {}
#
#        names_or_objs = {'M0': self.m0file, 'M1':self. m1file, 'M2': self.m2file, 'M3': self.m3file, 'M4': self.m4file}
#
#        # model components are interpolated to 2nd order
#        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
#            phase, wave, values = sncosmo.read_griddata_ascii(p2 + names_or_objs[key])
#            
#            values *= self._SCALE_FACTOR
#            # self._model[key] = Spline2d(phase, wave, values, kx=2, ky=2)
##            self._model[key] = values
#            self._model[key] = sncosmo.salt2utils.BicubicInterpolator(phase, wave, values)
#            # The "native" phases and wavelengths of the model are those
#            # of the first model component.
#            if key == 'M0':
#                self._phase = np.array([ -12.,  -9.,  -6.,  -3.,   0.,   3.,   6.,   9.,  12.,  15.,  18.,
#                                        21.,  24.,  27.,  30.,  33.,  36.,  39.,  42.,  45.,  48.])
#
#                                
#                self._wave = wave        
#            
#    def _flux(self, phase, wave):
#        """
#        """
#        m0 = self._model['M0'](phase, wave)
#        m1 = self._model['M1'](phase, wave)
#        m2 = self._model['M2'](phase, wave)
#        m3 = self._model['M3'](phase, wave)
#        m4 = self._model['M4'](phase, wave)
#        return (10. ** (-0.4 * (m0 + self._parameters[0] * m1 + self._parameters[1] * m2 + self._parameters[2] * m3 + self._parameters[3] * m4 + self._parameters[4] + 48.59)) / (wave ** 2 / 299792458. * 1.e-10))
#    
#
#        
#    def AstropyTable_flux(self, band_used = ['new_fI_10','fB_10','fV_10','fR_10','fU_10']):
#        """
#        """       
#        time = []
#        flux = []
#        fluxerr = []
#        band = []
#        zp = []
#        zpsys = []
#        vega = sncosmo.get_magsystem('vega_snf_10')
#        for i in range(len(self._phase)):
#            for b in band_used:
#            
#                filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+b+'.dat')
#                wlen = filt2[:,0]
#                tran = filt2[:,1]
#                self.splB = Spline1d(wlen, tran, k=1,ext = 1)    
#                
#                #computation of the integral
#                dt = 100000
#                xs = np.linspace(float(wl_min_sug), float(wl_max_sug), dt)
#                dxs = (float(wl_max_sug-wl_min_sug)/(dt-1))
#                spec_flux = self.model_spectrum_flux_m0(self._phase, xs)
##                plt.plot(xs,self.splB(xs))
##                plt.plot(xs,spec_flux[i,:])
##                print spec_flux[i,:] 
#                plt.show()
#                inte = np.sum(spec_flux[i,:]*(xs  / (CLIGHT*HPLANCK))*self.splB(xs)*dxs)
#                    
#                flux.append(inte)
#                fluxerr.append(0.00000000001)
#                time.append(self._phase[i])
#                band.append(b)
#                zp.append(2.5*np.log10(vega.zpbandflux(b)))
#                zpsys.append('vega_snf_10')
#        data = Table([time, band, flux, fluxerr, zp, zpsys], names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), meta={'name': 'data'})
#        return data