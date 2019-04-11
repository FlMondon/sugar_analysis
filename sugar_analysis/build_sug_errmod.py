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

from .constant import t_min_sug, t_max_sug 
from .fitting_lc import LC_Fitter
from .load_sugar import register_SUGAR
from .Hubble_fit import read_input_data_SNf 
from astropy.extern import six

def build_sig_par(nb_bin=9):
    bands = ['cspg', 'cspb', 'cspv', 'cspr', 'cspi']
    sigmas2 = []
    for band in bands:
        for j in range(nb_bin):
            sigmas2.append(band+'_'+str(j))
    
    return sigmas2
    
    
def make_method(obj):
    """Decorator to make the function a method of *obj*.
    In the current context::
      @make_method(Axes)
      def toto(ax, ...):
          ...
    makes *toto* a method of `Axes`, so that one can directly use::
      ax.toto()
    COPYRIGHT: from Yannick Copin
    """

    def decorate(f):
        setattr(obj, f.__name__, f)
        return f

    return decorate


def get_reml(output_path='../../sugar_analysis_data/err_mod_training/',
                 modeldir='../../sugar_model/', nb_bin=9, reml=False, fit_iter=True):
    """
    Parameters
    ----------
    """
    class build_sugar_error_model_case(build_sugar_error_model):
        freeparameters = build_sig_par(nb_bin=nb_bin)
    

    reml = build_sugar_error_model_case(output_path='../../sugar_analysis_data/err_mod_training/',
                 modeldir='../../sugar_model/', nb_bin=nb_bin, reml=reml, fit_iter=fit_iter)
    return reml

class build_sugar_error_model(object):
    """
    """
    
    def __new__(cls,*arg,**kwargs):
        """ Upgrade of the New function to enable the
        the _minuit_ black magic
        """
        obj = super(build_sugar_error_model,cls).__new__(cls)
        
        exec ("@make_method(build_sugar_error_model)\n"+\
             "def _minuit_reml_(self,%s): \n"%(", ".join(obj.freeparameters))+\
             "    sigmas2 = %s \n"%(", ".join(obj.freeparameters))+\
             "    return self.likelihood(sigmas2)\n")
        return obj
    
    
    def __init__(self, output_path='../../sugar_analysis_data/err_mod_training/',
                 modeldir='../../sugar_model/', nb_bin=9, reml=False, fit_iter=True):
        register_SUGAR(modeldir=modeldir, version='0.0')
        self.output_path = output_path
        self.source = sncosmo.get_source('sugar', version='0.0')
        self.nb_bin = nb_bin
        self.reml = reml
        self.fit_iter = fit_iter
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, 
                          effects=[dust], 
                          effect_names=['mw'], 
                          effect_frames=['obs']) 
        
        m0file = 'sugar_template_0.dat'
        m1file = 'sugar_template_1.dat'
        m2file = 'sugar_template_2.dat'
        m3file = 'sugar_template_3.dat'
        m4file = 'sugar_template_4.dat'

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
    
    def weight_matrix(self, sigmas2):
        
        band = self.dic[self.sn_name]['data_table']['band'][self.Filtre]
        data_fluxerr = self.dic[self.sn_name]['data_table']['fluxerr'][self.Filtre]
        data_flux = self.dic[self.sn_name]['data_table']['flux'][self.Filtre]
        cm_diag = np.zeros_like(data_fluxerr)
        phase_bin = np.linspace(t_min_sug, t_max_sug, self.nb_bin+1)
        
        
        
        for i in range(self.nb_bin):
            for j, b in enumerate(band):
                  phase_obs = self.dic[self.sn_name]['data_table']['time'][self.Filtre][j] - self.dic[self.sn_name]['res']['parameters'][1]
                  phase = phase_obs / (1 + self.dic[self.sn_name]['res']['parameters'][0])    
                  if b == 'cspg':
                     if phase >= phase_bin[i] and phase <= phase_bin[i+1]:
                         cm_diag[j] = 1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigmas2[i])
                         
                  elif b == 'cspb': 
                      if phase >= phase_bin[i] and phase <= phase_bin[i+1]:
                          cm_diag[j] =  1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigmas2[i+self.nb_bin])    
                          
                  elif b == 'cspv3014' or b == 'cspv9844' : 
                      if phase >= phase_bin[i] and phase <= phase_bin[i+1]:
                          cm_diag[j] =  1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigmas2[i+self.nb_bin*2])
                          
                  elif b == 'cspr': 
                      if phase >= phase_bin[i] and phase <= phase_bin[i+1]:
                          cm_diag[j] =  1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigmas2[i+self.nb_bin*3])         
                          
                  elif b == 'cspi': 
                      if phase >= phase_bin[i] and phase <= phase_bin[i+1]:
                          cm_diag[j] =  1/((data_fluxerr[j]*1.0857362047581294/data_flux[j] )**2 + sigmas2[i+self.nb_bin*4])
                  else:
                      raise ValueError ('filter have to be in this set [i, r, v3014, v9844, b, g]')
                  
        w = np.diag(cm_diag)
        det_cov = np.sum(np.log(cm_diag))
        return w, det_cov
    

    def likelihood(self, sigmas2):
        res = np.array([])
        w_m = []
        log_det_cov = 0.
        try:
            self.dic =  pkl.load(open(self.output_path+self.param_sug_path))
        except:
            self.dic = pkl.load(open(self.output_path+self.param_sug_path,
                                     'rb'), encoding='latin1')
            
        self.rids = read_input_data_SNf(res_dict_path=self.output_path+self.param_sug_path)
        self.dic = self.rids.delete_fit_fail(self.dic['data'])

        for sn_name in self.dic.keys():
            self.sn_name = sn_name
        

            self.Filtre = self.dic[self.sn_name]['res']['data_mask']    
            res = np.concatenate((res, self.residuals()), axis=None)
            w_i, log_det_cov_i = self.weight_matrix(sigmas2)
            w_m.append(w_i)
            log_det_cov += log_det_cov_i
        self.w = block_diag(*w_m)
        self.res = res
        if self.reml:
            self.H = self.build_H()
            counter_term = np.linalg.slogdet(np.dot(self.H.T,np.dot(self.w, self.H)))[1]
        else:
            counter_term = 0
        L = - log_det_cov + np.dot(self.res, np.dot(self.w, self.res)) + counter_term
        print L / len(self.res)
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
    
    def setup_guesses(self,**kwargs):
        """ Defines the guesses, boundaries and fixed values
        that will be passed to the given model.
        For each variable `v` of the model (see freeparameters)
        the following array will be defined and set to param_input:
           * v_guess,
           * v_boundaries,
           * v_fixed.
        Three arrays (self.paramguess, self.parambounds,self.paramfixed)
        will be accessible that will point to the defined array.
        Parameter
        ---------
        **kwargs the v_guess, v_boundaries and, v_fixed for as many
        `v` (from the freeparameter list).
        All the non-given `v` values will be filled either by pre-existing
        values in the model or with: 0 for _guess, False for _fixed, and
        [None,None] for _boundaries
        Return
        ------
        Void, defines param_input (and consquently paramguess, parambounds and paramfixed)
        """


        self.param_input = {}

        # -- Finally if no values have been set, let's do it
        for name in self.freeparameters:
                self.param_input[name+"_guess"] = 0.003
                self.param_input[name+"_boundaries"] = (0.,2.)
                
    def fit(self, **kwargs):
        """
        How to use kwargs 
        For each variable `v` of the model (see freeparameters)
        the following array will be defined and set to param_input:
           * v_guess,
           * v_boundaries,
           * v_fixed.
        Three arrays (self.paramguess, self.parambounds,self.paramfixed)
        will be accessible that will point to the defined array.
        Parameter
        ---------
        **kwargs the v_guess, v_boundaries and, v_fixed for as many
        `v` (from the freeparameter list).
        All the non-given `v` values will be filled either by pre-existing
        values in the model or with: 0 for _guess, False for _fixed, and
        [None,None] for _boundaries
        """
        #First step
        lcf = LC_Fitter(model_name='sugar', sample='csp',
                        modelcov=False, qual_crit=True, version_sug='0.0')
        lcf.fit_sample()
        self.param_sug_path = 'param_sugerrmod_0.pkl'
        lcf.write_result(specific_file=self.output_path+self.param_sug_path)
        self.setup_guesses(**kwargs)
        self._fit_minuit_()
        self.err_mod_path = 'train_intres_path_0.dat'
        self.write_res(self.err_mod_path)
        l = self._migrad_output_[0].fval / len(self.res)
        i = 0
        l_p = 0
        
        if self.fit_iter :
            while l < l_p and i <= 10 :
#                l_p = self._migrad_output_[0].fval / len(self.res)
                lcf = LC_Fitter(model_name='sugar', sample='csp',
                                modelcov=True, qual_crit=True, 
                                mod_errfile='../sugar_analysis_data/err_mod_training/'+self.err_mod_path, 
                                version_sug='%s.0'%str(i+1))
                lcf.fit_sample()
                self.param_sug_path = 'param_sugerrmod_%s.pkl'%str(i+1)
                lcf.write_result(specific_file=self.output_path+self.param_sug_path)
                self.setup_guesses(**kwargs)
                self._fit_minuit_()
                l = self._migrad_output_[0].fval / len(self.res)
                m_p = self._migrad_output_
                res_p = self.resultsfit
                i += 1
            self._migrad_output_ = m_p
            self.resultsfit = res_p
            self.err_mod_path = 'train_intres_path_%s.dat'%str(i+1)
            self.write_res(self.err_mod_path)
            
    def write_res(self, train_intres_path):
        
        train_err_mod = open(self.output_path+train_intres_path, 'w')
        csp_bands = ['cspb', 'cspg', 'cspv', 'cspr', 'cspi']
        for j, band in enumerate(csp_bands):
            if band == 'cspv':
                sncosmo_cspv9844 = sncosmo.get_bandpass('cspv9844')
                weff = sncosmo_cspv9844.wave_eff
                Warning('to complete cspv3014')
            else:
                sncosmo_band = sncosmo.get_bandpass(band)
                weff = sncosmo_band.wave_eff

            phase_bin = np.linspace(t_min_sug, t_max_sug, self.nb_bin+1)
            for i in range(self.nb_bin):
                    p_bin = np.min([phase_bin[i], phase_bin[i+1]])
                    value = self._migrad_output_[1][i+self.nb_bin*j]['value']
                    train_err_mod.write('%f %f %f \n'%(p_bin, weff, value))
                    
    
#    def intrinsic_dispertion(self, sigmas2):
#        return 
#    
                
    def _setup_minuit_(self):
        """
        """
        # == Minuit Keys == #
        minuit_kwargs = {}
        for param in self.freeparameters:
            minuit_kwargs[param] = self.param_input["%s_guess"%param]
            minuit_kwargs['limit_'+param] = self.param_input["%s_boundaries"%param]
        
        self.minuit =Minuit(self._minuit_reml_, **minuit_kwargs)
    
    def _fit_minuit_(self):
        """
        """
        self._setup_minuit_()
        self._migrad_output_ = self.minuit.migrad(ncall=10000)
        
        if self._migrad_output_[0]["is_valid"] is False:
            print("migrad is not valid")
            
            
        self.resultsfit = np.asarray([self.minuit.values[k]
                              for k in self.freeparameters])

                