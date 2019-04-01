#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 14:10:32 2019

@author: florian
"""

import sncosmo
import numpy as np
from scipy.linalg import block_diag

from iminuit import Minuit
try:
   import cPickle as pkl
except ModuleNotFoundError:
   import pickle as pkl
   
from .fitting_lc import LC_Fitter
from .load_sugar import register_SUGAR
from .Hubble_fit import read_input_data_SNf 

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
    
    def residuals(self):
        sys = sncosmo.get_magsystem('csp')
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
            data_mag[j] = sys.band_flux_to_mag(flux, band[j])
        model_mag = self.model.bandmag(band, 'csp', time_obs)
        
        residuals = data_mag - model_mag
        return residuals
    
    def cov_matrix(self, sigg, sigb, sigv, sigr, sigi):
        
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
                  raise ValueError ('filter have to be in this set [i, r, v3014, v9844, b, g ]')
        cov = np.diag(cm_diag)
        det_cov = np.sum(np.log(cm_diag))
        return cov, det_cov
    
    def restricted_likelihood(self, sigg, sigb, sigv, sigr, sigi):
        res = np.array([])
        cov_m = []
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
            cov_i, log_det_cov_i = self.cov_matrix(sigg, sigb, sigv, sigr, sigi)
            cov_m.append(cov_i)
            log_det_cov += log_det_cov_i
        cov_m = block_diag(*cov_m)
#        self.w = np.linalg.inv(cov_m)
        L = - log_det_cov + np.dot(res, np.dot(cov_m, res))
        return L
    
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