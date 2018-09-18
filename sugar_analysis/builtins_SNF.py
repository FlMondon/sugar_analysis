# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:34 2017
@author: mondon
Load SNF Filters and magnitude systeme and SUGAR
"""

import sncosmo 
import astropy.units as u
import numpy as np
from sncosmo.models import _SOURCES


sugar_model = '/home/florian/sugar_model/'
sugar_analysis_data = '/home/florian/sugar_analysis_data/'

def register_SNf_bands_width(width=10):
    
 
    band_file_U = 'USNf_3300-4102.dat'
    band_file_B = 'BSNf_4102-5100.dat'
    band_file_V = 'VSNf_5200-6289.dat'
    band_file_R = 'RSNf_6289-7607.dat'
    band_file_I = 'ISNf_7607-9200.dat'
       
        
    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/SNIFS/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/SNIFS/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/SNIFS/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/SNIFS/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/SNIFS/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]
    
    
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='USNf')    
    sncosmo.registry.register(band_U, force=True)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='BSNf')    
    sncosmo.registry.register(band_B, force=True)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='VSNf')    
    sncosmo.registry.register(band_V, force=True)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='RSNf')    
    sncosmo.registry.register(band_R, force=True)
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='ISNf')    
    sncosmo.registry.register(band_I, force=True)
    
    

    
    band_file_U = 'fU_'+str(width)+'.dat'
    band_file_B = 'fB_'+str(width)+'.dat'
    band_file_V = 'fV_'+str(width)+'.dat'
    band_file_R = 'fR_'+str(width)+'.dat'
    band_file_I = 'fI_'+str(width)+'.dat'
    band_file_new_U = 'new_fU_'+str(width)+'.dat'
    band_file_new_I = 'new_fI_'+str(width)+'.dat'
        
    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/Florian/'+ band_file_U)
    wl_U = filt2[:,0]
    transmission_U = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/Florian/'+ band_file_B)
    wl_B = filt2[:,0]
    transmission_B = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/Florian/'+ band_file_V)
    wl_V = filt2[:,0]
    transmission_V = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/Florian/'+ band_file_R)
    wl_R = filt2[:,0]
    transmission_R = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/Florian/'+ band_file_I)
    wl_I = filt2[:,0]
    transmission_I = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/Florian/'+ band_file_new_U)
    wl_new_U = filt2[:,0]
    transmission_new_U = filt2[:,1]

    filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/Florian/'+ band_file_new_I)
    wl_new_I = filt2[:,0]
    transmission_new_I = filt2[:,1]
            
    band_U = sncosmo.Bandpass(wl_U,transmission_U,wave_unit=u.AA,name='fU_'+str(width))    
    sncosmo.registry.register(band_U, force=True)
    band_B = sncosmo.Bandpass(wl_B,transmission_B,wave_unit=u.AA,name='fB_'+str(width))    
    sncosmo.registry.register(band_B, force=True)
    band_V = sncosmo.Bandpass(wl_V,transmission_V,wave_unit=u.AA,name='fV_'+str(width))    
    sncosmo.registry.register(band_V, force=True)
    band_R = sncosmo.Bandpass(wl_R,transmission_R,wave_unit=u.AA,name='fR_'+str(width))    
    sncosmo.registry.register(band_R, force=True)    
    band_I = sncosmo.Bandpass(wl_I,transmission_I,wave_unit=u.AA,name='fI_'+str(width))    
    sncosmo.registry.register(band_I, force=True)  
    band_new_U = sncosmo.Bandpass(wl_new_U,transmission_new_U, wave_unit=u.AA,name='new_fU_'+str(width))    
    sncosmo.registry.register(band_new_U, force=True) 
    band_new_I = sncosmo.Bandpass(wl_new_I,transmission_new_I, wave_unit=u.AA,name='new_fI_'+str(width))    
    sncosmo.registry.register(band_new_I, force=True)
    
def load_spectral_magsys_fits2(relpath, name=None, version=None):
    
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
    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0', load_spectral_magsys_fits2,
                             args=[sugar_analysis_data+'data/MagSys/bd_17d4708_stisnic_002.ascii'],
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
    sncosmo.registry.register_loader(sncosmo.MagSystem, 'vega_snf_0_'+str(width), load_spectral_magsys_fits2,
                             args=[sugar_analysis_data+'data/MagSys/bd_17d4708_stisnic_002.ascii'],
                             meta={'description': 'use bd_17d4708 spectrum that come from snfit'},
                             force=True)        

    bands_snf ={'fU_'+str(width): ('vega_snf_0_'+str(width), 9.787),
    'new_fU_'+str(width): ('vega_snf_0_'+str(width), 9.807),
    'new_fI_'+str(width): ('vega_snf_0_'+str(width),8.786 ),
	'fB_'+str(width): ('vega_snf_0_'+str(width), 9.791),
	'fV_'+str(width): ('vega_snf_0_'+str(width), 9.353),
	'fR_'+str(width): ('vega_snf_0_'+str(width), 9.011),
    'ISNf': ('vega_snf_0_'+str(width), 9.011),
    'fI_'+str(width): ('vega_snf_0_'+str(width), 8.768)}

    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_snf),'vega_snf_'+str(width), force=True) 

# =============================================================================
# Sources






class SUGARSource(sncosmo.Source):
    _param_names = ['Xgr', 'q1', 'q2', 'q3', 'A']
    param_names_latex = ['X_r', 'q_1', 'q_2', 'q_3', 'A']
    


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
        self._parameters = np.array([1.0e10-15, 1., 1., 1., 0.])
        




        names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}

        # model components are interpolated to 2nd order
        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
            phase, wave, values = sncosmo.read_griddata_ascii(sugar_model + names_or_objs[key])
            
            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
                
                self._phase = np.array([ -12.,  -9.,  -6.,  -3.,   0.,   3.,   6.,   9.,  12.,  15.,  18.,
                                        21.,  24.,  27.,  30.,  33.,  36.,  39.,  42.,  45.,  48.])
                                
                self._wave = wave
            self._model[key] = sncosmo.salt2utils.BicubicInterpolator(phase, wave, values) 


    def _flux(self, phase, wave):
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        m2 = self._model['M2'](phase, wave)
        m3 = self._model['M3'](phase, wave)
        m4 = self._model['M4'](phase, wave)
        return (self._parameters[0] * 10. ** (-0.4 * (m0 + self._parameters[1] * m1 + self._parameters[2] * m2 + self._parameters[3] * m3 + self._parameters[4] * m4 + 48.59)) / (wave ** 2 / 299792458. * 1.e-10))


    def bandflux_rcov(self, band, phase):
        return np.zeros(phase.shape, dtype=np.float64)
    
from sncosmo.builtins import DATADIR

# Sugar model
def load_sugarmodel(relpath, name=None, version=None):
    return SUGARSource(modeldir=relpath, name=name, version=version)

def register_SUGAR():
    from operator import itemgetter
    

    infile = open(sugar_model + 'SUGAR_model_v1.asci', 'r')
    lines = infile.readlines()
    infile.close()
    
    listik = []
    for line in lines:
        new_line = [float(i) for i in line.split()]
        listik.append(new_line)
    
    s = sorted(listik, key=itemgetter(0))
    
    names = ['0','1','2','3','4']
    for i in names:
        outfile = open(sugar_model + 'sugar_template_' + i + '.dat', 'w')
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
        _SOURCES.register_loader('sugar', load_sugarmodel, args=([sugar_model+'SUGAR_model_v1.asci']), version=ver, meta=meta, force=True)

   