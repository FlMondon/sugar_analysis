#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 14:10:32 2019

@author: florian
"""

import sncosmo
import numpy as np
import os
from scipy.linalg import block_diag

from iminuit import Minuit
try:
   import cPickle as pkl
except ModuleNotFoundError:
   import pickle as pkl
   
from .fitting_lc import LC_Fitter
from .load_sugar import register_SUGAR
from .Hubble_fit import read_input_data_SNf 
from astropy.extern import six


class build_sugar_error_model(object):
    """
    """
    def __init__(self, output_path='../../sugar_analysis_data',
                 modeldir='../../sugar_model/'):
        register_SUGAR(modeldir=modeldir)
        self.output_path = output_path
        self.source = sncosmo.get_source('sugar')
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, 
                          effects=[dust], 
                          effect_names=['mw'], 
                          effect_frames=['obs']) 
        
        m0file='sugar_template_0.dat'
        m1file='sugar_template_1.dat'
        m2file='sugar_template_2.dat'
        m3file='sugar_template_3.dat'
        m4file='sugar_template_4.dat'

        self.sys = sncosmo.get_magsystem('csp')
        names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}
        self.M_keys = ['M0', 'M1', 'M2', 'M3', 'M4']    
        # Make filenames into full paths.
        if modeldir is not None:
            for k in names_or_objs:
                v = names_or_objs[k]
                if (v is not None and isinstance(v, six.string_types)):
                    names_or_objs[k] = os.path.join(modeldir, v)
        self._model = {}        
        for i, key in enumerate(self.M_keys):
            phase, wave, values = sncosmo.read_griddata_ascii(names_or_objs[key])
            self._model[key] = sncosmo.salt2utils.BicubicInterpolator(phase, wave, values) 
    
    def residuals(self):
        
        self.model.set(z=float(self.dic[self.sn_name]['res']['parameters'][0]))
        self.model.set(mwebv=self.dic[self.sn_name]['mwebv'])
        self.model.set(q1=float(self.dic[self.sn_name]['res']['parameters'][3]))
        self.model.set(q2=float(self.dic[self.sn_name]['res']['parameters'][4]))
        self.model.set(q3=float(self.dic[self.sn_name]['res']['parameters'][5]))
        self.model.set(A=float(self.dic[self.sn_name]['res']['parameters'][6]))
        self.model.set(Xgr=float(self.dic[self.sn_name]['res']['parameters'][2]))
        self.model.set(t0=float(self.dic[self.sn_name]['res']['parameters'][1]))
        band = self.dic[self.sn_name]['data_table']['band'][self.Filtre]
        time_obs = self.dic[self.sn_name]['data_table']['time'][self.Filtre]
#        zp = self.dic[self.sn_name]['data_table']['zp'][self.Filtre]
        data_flux = self.dic[self.sn_name]['data_table']['flux'][self.Filtre]
        data_mag = np.zeros_like(data_flux)
        for j, flux in enumerate(data_flux):
            data_mag[j] = self.sys.band_flux_to_mag(flux, band[j])
        model_mag = self.model.bandmag(band, 'csp', time_obs)
        
        residuals = data_mag - model_mag
        return residuals
    
    def weight_matrix(self, sigg, sigb, sigv, sigr, sigi):
        
        band = self.dic[self.sn_name]['data_table']['band'][self.Filtre]
        data_fluxerr = self.dic[self.sn_name]['data_table']['fluxerr'][self.Filtre]
        data_flux = self.dic[self.sn_name]['data_table']['flux'][self.Filtre]
        cm_diag = np.zeros_like(data_fluxerr)
        for j, b in enumerate(band):
#              print b
              if b == 'cspg':
#                 print 'hello g'
                 cm_diag[j] = 1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigg)
              elif b == 'cspv3014' or b == 'cspv9844' : 
#                  print 'hello v'
                  cm_diag[j] =  1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigv)
              elif b == 'cspb': 
#                  print 'hello b'
                  cm_diag[j] =  1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigb)        
              elif b == 'cspr': 
#                  print 'hello r'
                  cm_diag[j] =  1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigr)         
              elif b == 'cspi': 
#                  print  'hello i'
                  cm_diag[j] =  1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigi)
              else:
                  raise ValueError ('filter have to be in this set [i, r, v3014, v9844, b, g]')
                  
        cov = np.diag(cm_diag)
        det_cov = np.sum(np.log(cm_diag))
        return cov, det_cov
    

    def restricted_likelihood(self, sigg, sigb, sigv, sigr, sigi):
        res = np.array([])
        w_m = []
        log_det_cov = 0.
        try:
            self.dic =  pkl.load(open(self.output_path+'param_sugerrmod.pkl'))
        except:
            self.dic = pkl.load(open(self.output_path+'param_sugerrmod.pkl',
                                     'rb'), encoding='latin1')
        self.rids = read_input_data_SNf(res_dict_path=self.output_path+'param_sugerrmod.pkl')
        self.dic = self.rids.delete_fit_fail(self.dic['data'])

        for sn_name in self.dic.keys():
            self.sn_name = sn_name
        

            self.Filtre = self.dic[self.sn_name]['res']['data_mask']    
            res = np.concatenate((res, self.residuals()), axis=None)
            w_i, log_det_cov_i = self.weight_matrix(sigg, sigb, sigv, sigr, sigi)
            w_m.append(w_i)
            log_det_cov += log_det_cov_i
        self.w = block_diag(*w_m)
        self.res = res
        self.H = self.build_H()
        counter_term = np.linalg.slogdet(np.dot(self.H.T,np.dot(self.w, self.H)))[1]
        L = - log_det_cov + np.dot(self.res, np.dot(self.w, self.res)) + counter_term
        return L

    def model_comp(self, band, phase, z):
        band = sncosmo.get_bandpass(band)
        wave = band.wave / (1+z)
        wave2_factor = (wave ** 2 / 299792458. * 1.e-10)
        comp_sug = {}
        for i, key in enumerate(self.M_keys):
                comp_sug[key] = np.sum((10. ** (-0.4 *self._model[key](phase, wave))/  wave2_factor)*band(band.wave)*np.min(band.dwave))
                comp_sug[key] = self.sys.band_flux_to_mag(comp_sug[key], band)
        return np.array([comp_sug['M0'], comp_sug['M1'], comp_sug['M2'], comp_sug['M3'], comp_sug['M4'], 1.])
        
    def build_H(self):
        H = np.ones((len(self.res), 6))
        n = 0 
        for sn_name in self.dic.keys():
            Filtre = self.dic[sn_name]['res']['data_mask']    
            for j, band in enumerate(self.dic[sn_name]['data_table']['band'][Filtre]):
                phase_obs = self.dic[sn_name]['data_table']['time'][Filtre][j] - self.dic[sn_name]['res']['parameters'][1]
                phase = phase_obs / (1 + self.dic[sn_name]['res']['parameters'][0])
                H[n] = self.model_comp(band, phase, self.dic[self.sn_name]['res']['parameters'][0])
                n +=1
        return H
    def compute_reml(self):
#        lcf = LC_Fitter(model_name='sugar', sample='csp',
#                        modelcov=False, qual_crit=True)
#        lcf.fit_sample()
#        lcf.write_result(specific_file=self.output_path+'param_sugerrmod.pkl')
        m = Minuit(self.restricted_likelihood, sigg=0.1,  sigb=0.1, sigv=0.1,
                   sigr=0.1, sigi=0.1, limit_sigg=(0., 3.), limit_sigb=(0., 3.),
                   limit_sigv=(0., 3.), limit_sigr=(0., 3.), limit_sigi=(0., 3.))
        m.migrad(ncall=10000)
        return m