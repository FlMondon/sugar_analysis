# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:34 2017

@author: mondon
"""

import sncosmo 
import astropy.units as u
import numpy as np
import pylab as P
from scipy.interpolate import SmoothBivariateSpline as Spline2d
from sncosmo.models import _SOURCES

from matplotlib import pyplot as plt
CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27        # [erg s]
p2 = '../sugar_model/'
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
 
    

def register_SNf_bands_width(width=10):
    
    band_file_U = 'fU_'+str(width)+'.dat'
    band_file_B = 'fB_'+str(width)+'.dat'
    band_file_V = 'fV_'+str(width)+'.dat'
    band_file_R = 'fR_'+str(width)+'.dat'
    band_file_I = 'fI_'+str(width)+'.dat'
    band_file_new_U = 'new_fU_'+str(width)+'.dat'
    band_file_new_I = 'new_fI_'+str(width)+'.dat'
        
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

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_new_U)
    wl_new_U = filt2[:,0]
    transmission_new_U = filt2[:,1]

    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/Florian/'+ band_file_new_I)
    wl_new_I = filt2[:,0]
    transmission_new_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_'+str(width))    
    sncosmo.registry.register(band_U)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_'+str(width))    
    sncosmo.registry.register(band_B)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_'+str(width))    
    sncosmo.registry.register(band_V)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_'+str(width))    
    sncosmo.registry.register(band_R)    
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='fI_'+str(width))    
    sncosmo.registry.register(band_I)  
    band_new_U = sncosmo.Bandpass(wl_new_U,transmission_new_U, wave_unit=u.AA,name='new_fU_'+str(width))    
    sncosmo.registry.register(band_new_U) 
    band_new_I = sncosmo.Bandpass(wl_new_I,transmission_new_I, wave_unit=u.AA,name='new_fI_'+str(width))    
    sncosmo.registry.register(band_new_I)
    
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


def mag_sys_SNF_width(width=10):
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
    'new_fU_'+str(width): ('vega_snf_0_'+str(width), 9.807),
    'new_fI_'+str(width): ('vega_snf_0_'+str(width),8.786 ),
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
    _SCALE_FACTOR = 1.


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
        self._parameters = np.array([1., 1., 1., 0., 37.])
        




        names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}

        # model components are interpolated to 2nd order
        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
            phase, wave, values = sncosmo.read_griddata_ascii(p2 + names_or_objs[key])
            
            values *= self._SCALE_FACTOR
            # self._model[key] = Spline2d(phase, wave, values, kx=2, ky=2)
            self._model[key] = values

            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
#                self._phase = np.array([ -12.,  -9.,  -6.,  -3.,   0.,   3.,   6.,   9.,  12.,  15.,  18.,
#                                        21.,  24.,  27.,  30.,  33.,  36.,  39.,  42.,  45.,  48.,  53., 57., 57.1, 57.2, 63., 67.])
                self._phase = np.array([ -12.,  -9.,  -6.,  -3.,   0.,   3.,   6.,   9.,  12.,  15.,  18.,
                                                        21.,  24.,  27.,  30.,  33.,  36.,  39.,  42.,  45.,  48., 48.5,  53., 57., 57.1, 57.2])
                                
                self._wave = wave

#    def _model_flux(self):
#        m0 = self._model['M0']
#        m1 = self._model['M1']
#        m2 = self._model['M2']
#        m3 = self._model['M3']
#        m4 = self._model['M4']
#        mod_flux = np.zeros([len(m0[:, 0])+6,len(m0[0])])
##        print mod_flux
#        for j in range(len(m0[0])):
#            
#            for i in range(len(m0[:, 0])):
#                mod_flux[i,j] = (10. ** (-0.4 * (m0[i,j] + self._parameters[0] * m1[i,j] + self._parameters[1] * m2[i,j] + self._parameters[2] * m3[i,j] + self._parameters[3] * m4[i,j] + self._parameters[4] + 48.59)) / (self._wave[j] ** 2 / 299792458. * 1.e-10))
##                mod_flux[i,j] = (10. ** (-0.4 * (m0[i,j] + self._parameters[0] * m1[i,j] + self._parameters[1] * m2[i,j] + self._parameters[2] * m3[i,j] + self._parameters[3] * m4[i,j] + self._parameters[4] + 48.59)) * (self._wave[j]  / (CLIGHT*HPLANCK)))
#        
#            mod_flux[len(m0[:, 0])+1,j] = mod_flux[len(m0[:, 0]),j]*0.75
#            mod_flux[len(m0[:, 0])+2,j] = mod_flux[len(m0[:, 0]),j]*0.25           
#            mod_flux[len(m0[:, 0])+3,j] = mod_flux[len(m0[:, 0]),j]*0.23
#            mod_flux[len(m0[:, 0])+4,j] = mod_flux[len(m0[:, 0]),j]*0.22
##        print (10. ** (-0.4 * (m0 + self._parameters[0] * m1 + self._parameters[1] * m2 + self._parameters[2] * m3 + self._parameters[3] * m4 + self._parameters[4] + 48.59)) / (self._wave ** 2 / 299792458. * 1.e-10))
#
#        return mod_flux
##                    
    def _model_flux(self):
        m0 = self._model['M0']
        m1 = self._model['M1']
        m2 = self._model['M2']
        m3 = self._model['M3']
        m4 = self._model['M4']
        mod_flux = np.zeros([len(m0[:, 0])+5,len(m0[0])])
#        print mod_flux
        for j in range(len(m0[0])):
            
            for i in range(len(m0[:, 0])):
                mod_flux[i,j] = (10. ** (-0.4 * (m0[i,j] + self._parameters[0] * m1[i,j] + self._parameters[1] * m2[i,j] + self._parameters[2] * m3[i,j] + self._parameters[3] * m4[i,j] + self._parameters[4] + 48.59)) / (self._wave[j] ** 2 / 299792458. * 1.e-10))
#                mod_flux[i,j] = (10. ** (-0.4 * (m0[i,j] + self._parameters[0] * m1[i,j] + self._parameters[1] * m2[i,j] + self._parameters[2] * m3[i,j] + self._parameters[3] * m4[i,j] + self._parameters[4] + 48.59)) * (self._wave[j]  / (CLIGHT*HPLANCK)))
            
            mod_flux[len(m0[:, 0])+1,j] = mod_flux[len(m0[:, 0]),j]*0.98
            mod_flux[len(m0[:, 0])+1,j] = mod_flux[len(m0[:, 0]),j]*0.75
            mod_flux[len(m0[:, 0])+2,j] = mod_flux[len(m0[:, 0]),j]*0.25           
            mod_flux[len(m0[:, 0])+3,j] = mod_flux[len(m0[:, 0]),j]*0.23
            mod_flux[len(m0[:, 0])+4,j] = mod_flux[len(m0[:, 0]),j]*0.22
        print (10. ** (-0.4 * (m0 + self._parameters[0] * m1 + self._parameters[1] * m2 + self._parameters[2] * m3 + self._parameters[3] * m4 + self._parameters[4] + 48.59)) / (self._wave ** 2 / 299792458. * 1.e-10))

        return mod_flux
                
    def _flux(self, phase, wave):



        flux_inter = sncosmo.salt2utils.BicubicInterpolator(self._phase, self._wave, self._model_flux())

        
        return flux_inter(phase,wave)



    def bandflux_rcov(self, band, phase):
        return np.zeros(phase.shape, dtype=np.float64)
    
from sncosmo.builtins import DATADIR

# Sugar model
def load_sugarmodel(relpath, name=None, version=None):
#    abspath = DATADIR.abspath(relpath, isdir=True)
    return SUGARSource(modeldir=relpath, name=name, version=version)

def register_SUGAR():
    from operator import itemgetter
    

    infile = open(p2 + 'SUGAR_model_v1.asci', 'r')
    lines = infile.readlines()
    infile.close()
    
    listik = []
    for line in lines:
        new_line = [float(i) for i in line.split()]
        listik.append(new_line)
    
    s = sorted(listik, key=itemgetter(0))
    
    names = ['0','1','2','3','4']
    for i in names:
        outfile = open(p2 + 'sugar_template_' + i + '.dat', 'w')
        for line in s:
            j = 2+int(i)
            outfile.write('%4.4f %8.8f %8.8f' %(line[0],line[1],line[j]))
            outfile.write('\n')
        outfile.close()
    website = 'http://no'
    PF16ref = ('PF16', 'PF et al. 2016 '
              '<http://arxiv.org/>')
    for topdir, ver, ref in [('SUGAR_model', '1.0', PF16ref)]:
        meta = {'type': 'SN Ia', 'subclass': '`~sncosmo.SUGARSource`',
                'url': website, 'reference': ref}
        _SOURCES.register_loader('sugar', load_sugarmodel, args=(['../sugar_model/SUGAR_model_v1.asci']), version=ver, meta=meta)

   