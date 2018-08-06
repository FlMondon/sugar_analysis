#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  3 19:36:21 2018

@author: florian
"""
import sugar_to_salt as sts
import numpy as np
p2 = '../sugar_model/'
p3 = '../sncosmo_jla/salt2-4/'
from operator import itemgetter
import sncosmo
import cPickle as pkl
from iminuit import Minuit as minuit
from sncosmo.salt2utils import BicubicInterpolator

class fit_t0_spectro(sts.sugar_spectrum):
    
    def __init__(self, modeldir=None,
                 m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
                    m2file='sugar_template_2.dat',
                    m3file='sugar_template_3.dat',
                 m4file='sugar_template_4.dat', 
                 name=None, version=None):
        

        self._SCALE_FACTOR = 1.
        self._param_names = ['q1', 'q2', 'q3', 'A', 'Mgr']
        
    
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
            
        
        self.step = 15
        self.name = name
        self.version = version
        self._model = {}
        self._parameters = np.array([0., 0., 0., 0., 37.])
        self._salt_params = [1., 0., 0.]
        names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}

        # model components are interpolated to 2nd order
        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
            phase, wave, values = sncosmo.read_griddata_ascii(p2 + names_or_objs[key])
            values *= self._SCALE_FACTOR
            # self._model[key] = Spline2d(phase, wave, values, kx=2, ky=2)
            self._model[key] = BicubicInterpolator(phase, wave, values)

            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
                self._phase = phase
                self._wave = wave
        self.dic_sts = None
        self.data_spectra = pkl.load(open('../sugar_analysis_data/spectra_snia.pkl'))
       
            
    def _model_flux(self,parameters, data_phase_rf, data_wave):
        """
        """
        self._parameters = parameters
        model_flux = self.model_spectrum_flux_m0(data_phase_rf, data_wave)
        list_model_flux = []
        for j in range(len(model_flux)):
            list_model_flux += list(model_flux[j,:])
            
        model_flux = np.array(list_model_flux)
        self._parameters = np.array([0., 0., 0., 0., 37.])
        return model_flux
    
    def flbda2fnu(self, x, y, var=None, backward=False):
        """Convert *x* [A], *y* [erg/s/cm2/A] to *y* [erg/s/cm2/Hz]. Se
        `var=var(y)` to get variance."""
        
        f = x**2 / 299792458. * 1.e-10 # Conversion factor                                                                                                                      
        if backward: 
            f = 1./f   
        if var is None:                # Return converted signa
            return y * f
        else:                          # Return converted variance
            return var * f**2
        
    def flbda2ABmag(self, x, y, ABmag0=48.59, var=None):
        """Convert *x* [A], *y* [erg/s/cm2/A] to `ABmag =                                                                                                                       
        -2.5*log10(erg/s/cm2/Hz) - ABmag0`. Set `var=var(y)` to get                                                                                                             
        variance."""          
        z = self.flbda2fnu(x,y)
        if var is None:
            return -2.5*np.log10(z) - ABmag0 
        else:
            return (2.5/np.log(10)/z)**2 * self.flbda2fnu(x,y,var=var)
                       
            
            
    def chi2(self, cst, q1, q2, q3, A, Mgr):
        """
        """
        sn_name = self.sn_name_t0_fit
        

        chi2 = 0
        for nspec in self.data_spectra[sn_name].keys():
            
            wave = self.data_spectra[sn_name][nspec]['X']
            parameters = np.array([q1, q2, q3, A, Mgr])
            
            phase = (self.data_spectra[sn_name][nspec]['phase_salt2'] + cst) 
            model_flux = self._model_flux(parameters, phase, wave)
            model_mag = self.flbda2ABmag(wave,model_flux) 
            data_flux = self.data_spectra[sn_name][nspec]['Y_flux_without_cosmology']
            data_mag = self.flbda2ABmag(wave,data_flux)
            err_mag = self.data_spectra[sn_name][nspec]['V']
            chi2 += np.sum((data_mag-model_mag)**2 / err_mag**2)
        print chi2
#        plt.plot(self.meta_spectra['SNF20070727-016']['spectra']['54330.2944444']['wave'],self.meta_spectra['SNF20070727-016']['spectra']['54330.2944444']['flux'])
#        plt.plot(self.meta_spectra['SNF20070727-016']['spectra']['54330.2944444']['wave'],model_flux, color ='r')
#        plt.show()

        return chi2

    def fit_spectro(self):
       """
       """
       self.sn_name_t0_fit = 'SNF20070727-016'
       m_object = minuit(self.chi2, cst=0.01, q1=0., q2=0., q3=0., A=0., Mgr=37.)
       m_object.migrad()
       return m_object
#        print a
        # Note that below we multiply by the scale factor to conserve
        # bolometric luminosity.
#        f = a * self._source._flux(phase, restwave)
