#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 13:02:35 2018

@author: florian
"""

import numpy as np
from operator import itemgetter
import sncosmo
from sugar_analysis import builtins as Build_SNF
from sncosmo.salt2utils import BicubicInterpolator
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
import copy
from astropy.table import Table
from sugar_analysis import math_toolbox as sam
import constant as cst
from matplotlib import pyplot as plt
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
import cPickle as pkl

Build_SNF.register_SNf_bands_width(width=10)
Build_SNF.mag_sys_SNF_width(width=10)
Build_SNF.register_SUGAR()
sugar_model = '/home/florian/sugar_model/'
sugar_analysis_data = '/home/florian/sugar_analysis_data/'
class sugar_simulation():
    
    def __init__(self, modeldir=None,
                 m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
                    m2file='sugar_template_2.dat',
                    m3file='sugar_template_3.dat',
                 m4file='sugar_template_4.dat', 
                 parameters_init=np.array([37., 0., 0., 0., 0.]),
                 name=None, version=None):
        

        self._SCALE_FACTOR = 1.
        self._param_names = ['q1', 'q2', 'q3', 'Av', 'grey']
        
    
        infile = open(sugar_model+ 'SUGAR_model_v1.asci', 'r')
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
            
        
        self.step = 15
        self.name = name
        self.version = version
        self._model_raw = {}
        self._model = {}
        self._parameters = parameters_init
        self._salt_params = [1., 0., 0.]
        names_or_objs = {'M0': m0file, 
                         'M1': m1file, 
                         'M2': m2file, 
                         'M3': m3file, 
                         'M4': m4file}
        
        # model components are interpolated to 2nd order
        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
            phase, wave, values = sncosmo.read_griddata_ascii(sugar_model 
                                                        + names_or_objs[key])
            values *= self._SCALE_FACTOR
            self._model_raw[key] = np.loadtxt(sugar_model+names_or_objs[key])
            self._model[key] = BicubicInterpolator(phase, wave, values)

            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
                self._phase = phase
                self._wave = wave
        self.dic_sts = None
        
        
    def model_spectrum_flux(self, phase, wave):
        '''
        '''
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        m2 = self._model['M2'](phase, wave)
        m3 = self._model['M3'](phase, wave)
        m4 = self._model['M4'](phase, wave)
#        return 10. ** (-0.4 * (m0 + 48.59)) / (wave** 2 / 299792458. * 1.e-10)
        return 10. ** (-0.4 * (m0 + self._parameters[1] * m1 + 
                               self._parameters[2] * m2 + 
                               self._parameters[3] * m3 +
                               self._parameters[4] * m4 +
                               self._parameters[0]  +
                               48.59))/ (wave** 2 / 299792458. * 1.e-10)


#        source_salt2 = sncosmo.get_source('salt2')
#        flux_salt2 = source_salt2._flux(phase, wave)
#        return flux_salt2  

    def spectrum_generator(self, parameters=None, phase=None, wave=None):
        """
        simulate sugar spectrum given a set of parameters
        parameters : array of sugar parameters (np.array([q1, q2, q3, A, Mgr]))
        mean_error : array of errors mean for different phase 
        """
        self.dic_spectrum = {}
        par_init = copy.deepcopy(self._parameters)
        if  not isinstance(phase, np.ndarray):
            phase = self._phase
        if  isinstance(parameters, np.ndarray):    
            self._parameters = parameters
        if  isinstance(wave, np.ndarray):
            wave = wave
        else:
            wave = np.linspace(cst.wl_min_sug,cst.wl_max_sug,197)
        flux = self.model_spectrum_flux(phase, wave)
        flux_err = np.zeros_like(flux)
        self.dic_spectrum['parameters'] = parameters
        self.dic_spectrum['wave'] = wave
        self.dic_spectrum['phase'] = phase
        self.dic_spectrum['fluxerr'] = flux_err 
        self.dic_spectrum['flux'] = flux            
        self._parameters = par_init 
        
    def error_generator_lc_factory(self, band):
        """
        """
        if band == 'new_fI_10' or band =='fI_10':
            error = np.abs(np.random.normal(0.0101, 0.0425, 1)[0])
        if band == 'fU_10':
            error = np.abs(np.random.normal(0.0146, 0.0741, 1)[0])
        if band == 'fB_10':
            error = np.abs(np.random.normal(0.0141, 0.06898, 1)[0])
        if band == 'fV_10':
            error = np.abs(np.random.normal(0.0139, 0.0646, 1)[0])
        if band == 'fR_10':
            error = np.abs(np.random.normal(0.0088, 0.0395, 1)[0])
        return error

    def integral_to_phot(self, flux, band):
        filt2 = np.genfromtxt(sugar_analysis_data+
                              'data/Instruments/Florian/'+band+'.dat')
        wlen = filt2[:,0]
        tran = filt2[:,1]
        self.splB = Spline1d(wlen, tran, k=1,ext = 1)    
        #computation of the integral
        dt = 10000
        dxs = (float(cst.wl_max_sug-cst.wl_min_sug)/(dt-1))
        inte = np.sum(flux*(self.xs/(cst.CLIGHT*cst.HPLANCK))*
                      self.splB(self.xs)*dxs)
        return inte
        
    def AstropyTable_flux(self, noise_level=None, phase=None,band_used=['new_fI_10','fB_10','fV_10','fR_10','fU_10']):
        """
        """    
        time = []
        flux = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('vega_snf_10')
        if phase == None:
            phase = self._phase
        for p in phase:
            for b in band_used:
                xs = np.linspace(float(cst.wl_min_sug), 
                                 float(cst.wl_max_sug), 10000)
                self.xs = xs            
                spec_flux = self.model_spectrum_flux(p, xs)
                phot_flux = self.integral_to_phot(spec_flux[0],b)  
                flux.append(phot_flux)
                time.append(p)
                band.append(b)
                zp.append(2.5*np.log10(vega.zpbandflux(b)))
                zpsys.append('vega_snf_10')
        if noise_level is None:
            flux_err = np.ones(len(flux)) * 1e-7
        else:
            noise = noise_level * max(flux)
            flux += np.random.normal(loc=0, scale=noise, size=len(flux))
            flux_err = np.ones(len(flux)) * noise
        data = Table([time, band, flux, flux_err, zp, zpsys], 
                     names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), 
                     meta={'name': 'data'})
        return data
                



    def get_param(self, sed, sed_err):
        
        W = np.eye(len(sed_err))*1./sed_err
        A = np.array([self._model_raw['M1'][:,2],
                      self._model_raw['M2'][:,2],
                      self._model_raw['M3'][:,2], 
                      self._model_raw['M4'][:,2], 
                      np.ones_like(self._model_raw['M0'][:,2])]).T
        left = np.linalg.inv(np.dot(A.T,np.dot(W, A)))
        right = np.dot(A.T, np.dot(W, sed - self._model_raw['M0'][:,2]))
        h = np.dot(left, right)
        return h
    
    def fit_spec_sugar(self):   
        spec_mag = []
        spec_mag_err = []
        for p in (self._phase):
            spec_mag_p = sam.flbda2ABmag(self._wave, 
                                         self.model_spectrum_flux(p, 
                                                                self._wave)[0])  
            spec_mag_p_err = np.ones_like(spec_mag_p)
            spec_mag += list(spec_mag_p)
            spec_mag_err += list(spec_mag_p_err)
        spec_mag = np.array(spec_mag)
        spec_mag_err = np.array(spec_mag_err)    
        h = self.get_param(spec_mag, spec_mag_err)
        return h
    
        
    def sample_genarator(self, nb_gen):
        '''
        '''
        params = pkl.load(open(sugar_model+'sugar_parameters.pkl'))
        coefs = []
        for param_name in self._param_names:
            coefs.append([v[param_name] for v in params.values()])
        coefs = np.array(coefs).T
        coefs.shape
        self.coefs = coefs
        kde = KernelDensity(bandwidth=0.1)
        kde.fit(self.coefs)       
        v, l, vinv = np.linalg.svd(np.cov(coefs.T))
        rescaled_coefs = (np.dot(np.dot(np.diag(1./np.sqrt(l)),
                                        vinv), 
                                        coefs.T)).T    
        params = {'bandwidth': np.logspace(-2, 1, 300)}
        grid = GridSearchCV(KernelDensity(), params)
        grid.fit(rescaled_coefs)
        kde = grid.best_estimator_
        samples = kde.sample(nb_gen)
        kde_info = {'rot': v,
                    'rot_inv': vinv,
                    'eigenvals': l,
                    'kde': kde}
        self.samples = (np.dot(np.dot(kde_info['rot'], 
                                      np.diag(np.sqrt(kde_info['eigenvals']))), 
                                      samples.T)).T      
        
                
