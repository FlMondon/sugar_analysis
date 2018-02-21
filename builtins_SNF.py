# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:34 2017

@author: mondon
"""

import sncosmo 
import astropy.units as u
import numpy as np
import pylab as P
from scipy.interpolate import RectBivariateSpline as Spline2d
from sncosmo.models import _SOURCES

CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27        # [erg s]

def register_SNf_bands():
    
#    wl = np.arange(3000,9500,0.1)
#    transmission_U = np.zeros(len(wl))
#    transmission_B = np.zeros(len(wl))
#    transmission_V = np.zeros(len(wl))
#    transmission_R = np.zeros(len(wl))
#    transmission_I = np.zeros(len(wl))

#    for i in range (len(wl)): 
#        if wl[i] >= 3300 and wl[i] <= 4102:
#            transmission_U[i]= 1.
#        if wl[i] >= 4102 and wl[i] <= 5100:
#            transmission_B[i]= 1. 
#        if wl[i] >= 5200 and wl[i] <= 6289:
#            transmission_V[i]= 1.
#        if wl[i] >= 6289 and wl[i] <= 7607 :
#            transmission_R[i]= 1.
#        if wl[i] >= 7607 and wl[i] <= 9200 :
#            transmission_I[i]= 1.        
#    band_U = sncosmo.Bandpass(wl,transmission_U,wave_unit=u.AA,name='USNf')    
#    sncosmo.registry.register(band_U)
#    band_B = sncosmo.Bandpass(wl,transmission_B,wave_unit=u.AA,name='BSNf')    
#    sncosmo.registry.register(band_B)
#    band_V = sncosmo.Bandpass(wl,transmission_V,wave_unit=u.AA,name='VSNf')    
#    sncosmo.registry.register(band_V)
#    band_R = sncosmo.Bandpass(wl,transmission_R,wave_unit=u.AA,name='RSNf')    
#    sncosmo.registry.register(band_R)
#    band_I = sncosmo.Bandpass(wl,transmission_I,wave_unit=u.AA,name='ISNf')    
#    sncosmo.registry.register(band_I)
            
    band_file_U = 'USNf_3300-4102.dat'
    band_file_B = 'BSNf_4102-5100.dat'
    band_file_V = 'VSNf_5200-6289.dat'
    band_file_R = 'RSNf_6289-7607.dat'
    band_file_I = 'ISNf_7607-9200.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
    
    
    


#            
#    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='USNf',top_hat=True)    
#    sncosmo.registry.register(band_U)
#    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='BSNf',top_hat=True)    
#    sncosmo.registry.register(band_B)
#    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='VSNf',top_hat=True)    
#    sncosmo.registry.register(band_V)
#    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='RSNf',top_hat=True)    
#    sncosmo.registry.register(band_R)
#    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='ISNf',top_hat=True)    
#    sncosmo.registry.register(band_I)
#    
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='USNf')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='BSNf')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='VSNf')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='RSNf')    
    sncosmo.registry.register(band_R)
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='ISNf')    
    sncosmo.registry.register(band_I)
    
    
    band_file_U = 'fU.dat'
    band_file_B = 'fB.dat'
    band_file_V = 'fV.dat'
    band_file_R = 'fR.dat'

       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]


            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU')    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB')    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV')    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR')    
    sncosmo.registry.register(band_R)
 
    

def register_SNf_bands_width(width=20):
    
    band_file_U = 'fU_'+str(width)+'.dat'
    band_file_B = 'fB_'+str(width)+'.dat'
    band_file_V = 'fV_'+str(width)+'.dat'
    band_file_R = 'fR_'+str(width)+'.dat'
    band_file_I = 'fI_'+str(width)+'.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_'+str(width))    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_'+str(width))    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_'+str(width))    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_'+str(width))    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_R,transmission_I,wave_unit=u.AA,name='fI_'+str(width))    
    sncosmo.registry.register(band_I)  
            
def load_spectral_magsys_fits2(relpath, name=None, version=None):
#    fits = pyfits.open(relpath)
#    fit = fits[1]
#    wl = np.zeros(len(fit.data))
#    flux = np.zeros(len(fit.data))
#    for i in range(len(fit.data)):
#         wl[i] = fit.data[i][0]
#         flux[i] = fit.data[i][1]
#    refspectrum = sncosmo.spectrum.Spectrum(wl, flux,
#                           unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    
    data = np.genfromtxt(relpath) 
    dispersion = data[:,0]
    flux_density = data[:,1]
    refspectrum = sncosmo.spectrum.Spectrum(dispersion, flux_density,
                           unit=(u.erg / u.s / u.cm**2 / u.AA), wave_unit=u.AA)
    return sncosmo.magsystems.SpectralMagSystem(refspectrum, name=name)
    

def mag_sys_SNF():
    """
    define magnitude systeme for snf
    """
    
#    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0', load_spectral_magsys_fits2,
#                             args=['../sugar_analysis_data/Vega.fits'],
#                             meta={'description': 'use vega spectrum that come from snf data'}) 
#    bands_snf ={'USNF': ('vega_snf_0', 0.),
#	'BSNF': ('vega_snf_0', 0.),
#	'VSNF': ('vega_snf_0', 0.),
#	'RSNF': ('vega_snf_0', 0.),
#	'ISNF': ('vega_snf_0', 0)}
    
    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0', load_spectral_magsys_fits2,
                             args=['../sugar_analysis_data/data/MagSys/bd_17d4708_stisnic_002.ascii'],
                             meta={'description': 'use bd_17d4708 spectrum that come from snfit'})        

    bands_snf ={'USNf': ('vega_snf_0', 9.787),
	'BSNf': ('vega_snf_0', 9.791),
	'VSNf': ('vega_snf_0', 9.353),
	'RSNf': ('vega_snf_0', 9.011),
	'ISNf': ('vega_snf_0', 8.768),
    'fU': ('vega_snf_0', 9.787),
	'fB': ('vega_snf_0', 9.791),
	'fV': ('vega_snf_0', 9.353),
	'fR': ('vega_snf_0', 9.011)}

    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_snf),'vega_snf') 


def mag_sys_SNF_width(width=20):
    """
    define magnitude systeme for snf
    """
    
#    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0', load_spectral_magsys_fits2,
#                             args=['../sugar_analysis_data/Vega.fits'],
#                             meta={'description': 'use vega spectrum that come from snf data'}) 
#    bands_snf ={'USNF': ('vega_snf_0', 0.),
#	'BSNF': ('vega_snf_0', 0.),
#	'VSNF': ('vega_snf_0', 0.),
#	'RSNF': ('vega_snf_0', 0.),
#	'ISNF': ('vega_snf_0', 0)}
    
    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0_'+str(width), load_spectral_magsys_fits2,
                             args=['../sugar_analysis_data/data/MagSys/bd_17d4708_stisnic_002.ascii'],
                             meta={'description': 'use bd_17d4708 spectrum that come from snfit'})        

    bands_snf ={'fU_'+str(width): ('vega_snf_0_'+str(width), 9.787),
	'fB_'+str(width): ('vega_snf_0_'+str(width), 9.791),
	'fV_'+str(width): ('vega_snf_0_'+str(width), 9.353),
	'fR_'+str(width): ('vega_snf_0_'+str(width), 9.011),
    'fI_'+str(width): ('vega_snf_0_'+str(width), 8.768)}

    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_snf),'vega_snf_'+str(width)) 

 
def plot_SNf_filters(wl, transmission_U, transmission_V, transmission_B, transmission_R, transmission_I):
    P.plot(wl, transmission_U)
    P.plot(wl, transmission_B)
    P.plot(wl, transmission_V)
    P.plot(wl, transmission_R)
    P.plot(wl, transmission_I)
    P.ylim(0,2)
    P.show()

# =============================================================================
# Sources

class SUGARSource(sncosmo.Source):

    _param_names = ['q1', 'q2', 'q3', 'A', 'Mgr']
    param_names_latex = ['q_1', 'q_2', 'q_3', 'A', 'M_g']
    _SCALE_FACTOR = 1


    def __init__(self, modeldir=None,
                 m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
		 m2file='sugar_template_2.dat',
		 m3file='sugar_template_3.dat',
		 m4file='sugar_template_4.dat',
                 name=None, version=None):
        self.name = name
        self.version = version
        self._model = {}
        self._parameters = np.array([1., 1., 1., 1., 40.])

        print modeldir 	
        phase, wave, M0, M1, M2, M3 , M4, F = np.genfromtxt(modeldir) 
        model_ascii = {'M0': M0, 'M1': M1, 'M2': M2, 'M3': M3, 'M4': M4}
        

        # model components are interpolated to 2nd order
        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
            values = model_ascii[key]
            values *= self._SCALE_FACTOR
            self._model[key] = Spline2d(phase, wave, values, kx=2, ky=2)

            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
                self._phase = phase
                self._wave = wave


    def _flux(self, phase, wave):
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        m2 = self._model['M2'](phase, wave)
        m3 = self._model['M3'](phase, wave)
        m4 = self._model['M4'](phase, wave)
	return (10. ** (-0.4 * (m0 + self._parameters[0] * m1 + self._parameters[1] * m2 + self._parameters[2] * m3 + self._parameters[3] * m4 + self._parameters[4] + 48.59)) / (wave ** 2 / 299792458. * 1.e-10))
	#return (10. ** (-0.4 * (m0 + self._parameters[0] * m1  + 48.59)) / (wave ** 2 / 299792458. * 1.e-10))/ 10**self._parameters[4]
	#return (10. ** (-0.4 * (m0 + self._parameters[0] * m1 + self._parameters[1] * m2 + self._parameters[2] * m3 + self._parameters[3] * m4 + self._parameters[4]))/f)

        #Flux_nu=10**(-0.4*(self.Y_new_binning+ABmag0))
        #f = self.SUGAR_Wavelength**2 / 299792458. * 1.e-10
        #self.Flux=Flux_nu/f



from sncosmo.builtins import DATADIR

# Sugar model
def load_sugarmodel(relpath, name=None, version=None):
#    abspath = DATADIR.abspath(relpath, isdir=True)
    return SUGARSource(modeldir=relpath, name=name, version=version)

def register_SUGAR():
    website = 'http://no'
    PF16ref = ('PF16', 'PF et al. 2016 '
              '<http://arxiv.org/>')
    for topdir, ver, ref in [('SUGAR_model', '1.0', PF16ref)]:
        meta = {'type': 'SN Ia', 'subclass': '`~sncosmo.SUGARSource`',
                'url': website, 'reference': ref}
        _SOURCES.register_loader('sugar', load_sugarmodel, args=(['../sugar_analysis_data/SUGAR_model_v1.asci']), version=ver, meta=meta)

   