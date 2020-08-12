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
from scipy.interpolate import interp2d

from iminuit import Minuit
try:
   import cPickle as pkl
except ModuleNotFoundError:
   import pickle as pkl

from .constant import t_min_sug, t_max_sug 
from .fitting_lc import LC_Fitter
from .load_sugar import register_SUGAR
from .load_salt2_newerr import register_salt2_newerr
from .Hubble_fit import read_input_data_SNf 
from astropy.extern import six

try:
    import george 
    from george import kernels
    from scipy.optimize import minimize
    
except ModuleNotFoundError:
    Warning('george needed to run build error sugar gp')
def build_sig_par(nb_node=9):
    bands = ['cspb', 'cspg', 'cspv', 'cspr', 'cspi']
    sigmas2 = []
    for band in bands:
        for j in range(nb_node):
            sigmas2.append(band+'_'+str(j))
    
    return sigmas2, band
    
    
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


def get_reml(res_dict_path="../../sugar_analysis_data/resfitlc_csp_drop['cspu']_sugar_gold.pkl",
             output_path='../../sugar_analysis_data/err_mod_training/',
                 modeldir='../../sugar_model/', nb_node=4, sad_path = '../../',
                 fit_spline=True, reml=False, fit_iter=True, grid=None, model='sugar', select=None):
    """
    Parameters
    ----------
    """
    
    class build_sugar_error_model_case(build_sugar_error_model):
        freeparameters, bands = build_sig_par(nb_node=nb_node)
    

    reml = build_sugar_error_model_case(res_dict_path=res_dict_path, output_path=output_path, 
                 modeldir=modeldir, nb_node=nb_node, sad_path=sad_path, 
                 reml=reml, fit_iter=fit_iter, fit_spline=True, grid=grid, model=model, select=select)
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
    
    
    def __init__(self, res_dict_path="../../sugar_analysis_data/resfitlc_csp_drop['cspu']_sugar_gold.pkl",
                 output_path='../../sugar_analysis_data/err_mod_training/',
                 modeldir='../../sugar_model/', nb_node=4, reml=False, sad_path = '../../',
                 fit_iter=True, bands=['cspb', 'cspg', 'cspv', 'cspr', 'cspi'],
                 fit_spline=True, grid=None, model='sugar', select=None):
        
        self.modeldir = modeldir
        self.sad_path = sad_path
        if grid == None:
            self.phase_bin = np.linspace(t_min_sug, t_max_sug, self.nb_node)
        else:
            if len(grid) != nb_node:
                raise ValueError('len(grid) have to be equal to nb_node ')
            else:
                self.phase_bin = grid
        self.output_path = output_path
        self.bands = bands
        dic = pkl.load(open(res_dict_path,'rb'), encoding='latin1')
        self.sys = sncosmo.get_magsystem('csp')
        self.nb_node = nb_node
        self.nb_node = nb_node
        self.reml = reml
#        self.val_sample = os.listdir(sad_path+'sugar_analysis_data/DR3/validation')
#        self.training_sample = os.listdir(sad_path+'sugar_analysis_data/DR3/trainning_all')
#        Warning('training is all sample')
        dic = pkl.load(open(res_dict_path,'rb'))
        dic_res = dic['data']
        if type(select) != type(None):
            dic = {}
            for sn_name in select:
                dic[sn_name] = dic_res[sn_name]
            if len(dic) != 0:
                self.dic = dic
            else:
                raise ValueError('select not compatible with dic res')
        else:
            self.dic = dic_res
        self.like_val_it = [] #values of likelihood/N after each iteration
        self.fit_spline = fit_spline
        self.fit_iter = fit_iter
        self.v = 0.0
        print(model)
        self.model_name = model
        register_SUGAR(modeldir=self.modeldir, mod_errfile=None)
        if self.model_name == 'sugar':
            self.source = sncosmo.get_source('sugar')
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
    
            
            names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}
            self.M_keys = ['M0', 'M1', 'M2', 'M3', 'M4']    
            # Make filenames into full paths.
            if modeldir is not None:
                for k in names_or_objs:
                    v = names_or_objs[k]
                    names_or_objs[k] = os.path.join(modeldir, v)
            self._model = {}        
            for i, key in enumerate(self.M_keys):
                phase, wave, values = sncosmo.read_griddata_ascii(names_or_objs[key])
                self._model[key] = sncosmo.salt2utils.BicubicInterpolator(phase, wave, values) 
        else:
            register_salt2_newerr(modeldir=modeldir, mod_errfile=modeldir+'err_mod_4node.dat',
                       version='0.0')
            self.source = sncosmo.get_source('salt2_newerr', version='0.0')
            dust = sncosmo.CCM89Dust()
            self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs'])             

    def change_model_params(self, sn_name):
        
        self.model.set(z=float(self.dic[sn_name]['res']['parameters'][0]))
        self.model.set(mwebv=self.dic[sn_name]['mwebv'])
        self.model.set(t0=float(self.dic[sn_name]['res']['parameters'][1]))
        if self.model_name == 'sugar' :
            self.model.set(q1=float(self.dic[sn_name]['res']['parameters'][3]))
            self.model.set(q2=float(self.dic[sn_name]['res']['parameters'][4]))
            self.model.set(q3=float(self.dic[sn_name]['res']['parameters'][5]))
            self.model.set(A=float(self.dic[sn_name]['res']['parameters'][6]))
            self.model.set(Xgr=float(self.dic[sn_name]['res']['parameters'][2]))
        elif self.model_name == 'salt2_newerr':
            self.model.set(x1=float(self.dic[sn_name]['res']['parameters'][3]))
            self.model.set(c=float(self.dic[sn_name]['res']['parameters'][4]))
            self.model.set(x0=float(self.dic[sn_name]['res']['parameters'][2]))      
            
    def residuals(self, sn_name):   
        self.change_model_params(sn_name)
        Filtre = self.dic[sn_name]['res']['data_mask']         
        band = self.dic[sn_name]['data_table']['band'][Filtre]
        time_obs = self.dic[sn_name]['data_table']['time'][Filtre]
        phase = (time_obs- self.dic[sn_name]['res']['parameters'][1]) / (1 + self.dic[sn_name]['res']['parameters'][0])  
        data_flux = self.dic[sn_name]['data_table']['flux'][Filtre]
        model_mag = self.model.bandmag(band, 'csp', time_obs)
        mag = np.zeros_like(model_mag)
        for i, d in enumerate(data_flux):
            mag[i] = self.sys.band_flux_to_mag(d, band[i])

#        residuals = data_flux - model_flux
        residuals = (mag - model_mag)
        return residuals,  band, phase
    
    def weight_matrix(self, sn_name):
        Filtre = self.dic[sn_name]['res']['data_mask'] 
        
        band = self.dic[sn_name]['data_table']['band'][Filtre]
        data_fluxerr = self.dic[sn_name]['data_table']['fluxerr'][Filtre]
        data_flux = self.dic[sn_name]['data_table']['flux'][Filtre]
        self.change_model_params(sn_name)
        cm_diag = np.zeros_like(data_fluxerr)
        c = np.zeros_like(data_fluxerr)
        for j, b in enumerate(band):
              phase_obs = self.dic[sn_name]['data_table']['time'][Filtre][j] - self.dic[sn_name]['res']['parameters'][1]
              phase = phase_obs / (1 + self.dic[sn_name]['res']['parameters'][0])   
              if b == 'cspv3014' or b == 'cspv9844':
                    f = sncosmo.get_bandpass('cspv9844')
              else :
                    f = sncosmo.get_bandpass(b)
                    
              weff = f.wave_eff/(1 + self.dic[sn_name]['res']['parameters'][0])    
              data_magerr = (data_fluxerr * 1.0857362047581294) / data_flux
              c[j] = ((data_magerr[j])**2 + (self.intrinsic_dispertion(phase, weff)/ 1.0857362047581294)**2) #flux ???
              cm_diag[j] = 1/c[j]
                  
        w = cm_diag
        det_cov = np.sum(np.log(c))
        return w, det_cov
                    
                    
    
    def intrinsic_dispertion(self, p, wl):
        return self.spline(p, wl)[0]
    
    def likelihood(self, sigmas2):
        chi2 = 0.
        log_det_cov = 0.
        self.wave_bin = np.zeros(len(self.bands))
        node_array =  np.zeros((len(self.wave_bin),len(self.phase_bin)))
        for l in range(len(self.phase_bin)):   
            for i in range(len(self.wave_bin)):
                if l ==0:
                    if self.bands[i] == 'cspv':
                        f = sncosmo.get_bandpass('cspv9844')
                    else :
                        f = sncosmo.get_bandpass(self.bands[i])
                    self.wave_bin[i] = float(f.wave_eff)
                    
                node_array[i,l] = sigmas2[l+len(self.phase_bin)*i]
        self.spline = interp2d(self.phase_bin, self.wave_bin, node_array) 
        if self.reml:
            w_m = []
            nb_p = 0
        for sn_name in self.dic.keys():  
            if self.fit_spline:
                w_i, log_det_cov_i = self.weight_matrix(sn_name)
            else:
                w_i, log_det_cov_i = self.weight_matrix_bin(sn_name)
            log_det_cov += log_det_cov_i
            chi2 += np.sum(self.res_dic[sn_name]['res']**2*w_i)
            if self.reml:
                w_m.append(np.diag(w_i))
                nb_p += len(w_i)
        if self.reml:
            self.w = block_diag(*w_m)
            self.H = self.build_H(nb_p)
            counter_term = np.linalg.slogdet(np.dot(self.H.T,np.dot(self.w, self.H)))[1]
        else:
            counter_term = 0
        L = - log_det_cov + chi2 + counter_term
        print(-L, log_det_cov ,chi2 , counter_term)
        return L

    def model_comp(self, band, phase, z):
        band = sncosmo.get_bandpass(band)
        wave, dwave = sncosmo.utils.integration_grid(band.minwave(), band.maxwave(),
                                       5.0)
        wave = wave / (1+z)
        wave2_factor = (wave ** 2 / 299792458. * 1.e-10)
        comp_sug = {}
        for i, key in enumerate(self.M_keys):
                comp_sug[key] = np.sum((10. ** (-0.4 *self._model[key](phase, wave))/  wave2_factor)*band(wave)*dwave)
                comp_sug[key] = self.sys.band_flux_to_mag(comp_sug[key], band)
        return np.array([comp_sug['M0'], comp_sug['M1'], comp_sug['M2'], comp_sug['M3'], comp_sug['M4'], 1.])
        
    def build_H(self, nb_p):
        H = np.ones((nb_p, 6))
        n = 0 
        for sn_name in self.dic.keys():
            Filtre = self.dic[sn_name]['res']['data_mask']    
            for j, band in enumerate(self.dic[sn_name]['data_table']['band'][Filtre]):
                phase_obs = self.dic[sn_name]['data_table']['time'][Filtre][j] - self.dic[sn_name]['res']['parameters'][1]
                phase = phase_obs / (1 + self.dic[sn_name]['res']['parameters'][0])
                H[n] = self.model_comp(band, phase, self.dic[sn_name]['res']['parameters'][0])
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
                self.param_input[name+"_guess"] = 0.07
                self.param_input[name+"_boundaries"] = (0.,2.)
                
    def fit(self, csp_file="../../sugar_analysis_data/resfitlc_csp_drop['cspu']_sugar_gold.pkl", **kwargs):
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
        self.nb_point = 0.
        self.res_dic = {}
        for sn_name in self.dic.keys():
            self.res_dic[sn_name] = {}
            self.res_dic[sn_name]['res'], self.res_dic[sn_name]['band'], self.res_dic[sn_name]['p'] = self.residuals(sn_name)
            self.nb_point += len(self.res_dic[sn_name])
        self.setup_guesses(**kwargs)
        self._fit_minuit_()
        self.err_mod_path = 'train_intres_0_%snode.dat'%str(self.nb_node)
        self.write_res(self.err_mod_path)
        l = self._migrad_output_[0].fval 
        self.like_val_it.append(l)
        l_p = 0
        i = 0

        if self.fit_iter :
            while l < l_p and i <= 10 :
                l_p = l
                lcf = LC_Fitter(model_name=self.model_name, sample='csp', sad_path=self.sad_path,
                                modelcov=True, qual_crit=False, 
                                mod_errfile=self.sad_path+'sugar_analysis_data/err_mod_training/'+self.model_name+'/'+self.err_mod_path, 
                                version_sug=str(self.v+2),  modeldir = self.modeldir, sub_sample=list(self.dic.keys()), 
                                filter_drop_csp = ['cspu'], t0_fix=True, csp_file=csp_file)
                lcf.fit_sample()
                self.param_sug_path = 'param_errmod_%sit_%snode.pkl'%(str(i+1), str(self.nb_node))
                lcf.write_result(specific_file=self.output_path+self.model_name+'/'+self.param_sug_path)
                self.res_dic = {}
                try:
                    self.dic =  pkl.load(open(self.output_path+self.model_name+'/'+self.param_sug_path))
                except:
                    self.dic = pkl.load(open(self.output_path+self.model_name+'/'+self.param_sug_path,
                                             'rb'), encoding='latin1')
                self.dic = self.dic['data']
                self.nb_point = 0.
                for sn_name in self.dic.keys():
                    self.res_dic[sn_name] = {}
                    self.res_dic[sn_name]['res'], self.res_dic[sn_name]['band'], self.res_dic[sn_name]['p'] = self.residuals(sn_name)
                self.setup_guesses(**kwargs)
                self._fit_minuit_()
                l = self._migrad_output_[0].fval 
                self.like_val_it.append(l)
                m_p = self._migrad_output_
                res_p = self.resultsfit
                self.v += 1
                i +=1
                self.err_mod_path = 'train_intres_%s_%snode.dat'%(str(i), str(self.nb_node))
                self.write_res(self.err_mod_path)
            self._migrad_output_ = m_p
            self.resultsfit = res_p
            self.err_mod_path = 'err_mod_%snode.dat'%str(self.nb_node)
            self.write_res(self.err_mod_path)
            print(self.nb_point)
    def write_res(self, train_intres_path):
        
        train_err_mod = open(self.output_path+self.model_name+'/'+train_intres_path, 'w')
        if self.fit_spline:
            for i, p_bin in enumerate(self.phase_bin):
                for j, weff in enumerate(self.wave_bin):
                    value = self._migrad_output_[1][i+self.nb_node*j]['value']
                    train_err_mod.write('%f %f %f \n'%(p_bin, weff, value))   
        else:
            for i in range(self.nb_node):
                p_bin = np.mean([self.phase_bin[i], self.phase_bin[i+1]])
                for j, band in enumerate(self.bands):
                    if band == 'cspv':
                        sncosmo_cspv9844 = sncosmo.get_bandpass('cspv9844')
                        weff = sncosmo_cspv9844.wave_eff
                        Warning('to complete cspv3014')
                    else:
                        sncosmo_band = sncosmo.get_bandpass(band)
                        weff = sncosmo_band.wave_eff
                    value = self._migrad_output_[1][i+self.nb_node*j]['value']
                    train_err_mod.write('%f %f %f \n'%(p_bin, weff, value))
                
    def best_nb_node(self, err_path, csp_file=None, val_sample=None):
        
        best_node = None
        best_likelihood = 0.
        list_err_mod = os.listdir(err_path)
        self.v = 5
        for err_mod in list_err_mod:
            if err_mod.startswith('err_mod_'):
                self.v +=1
                err_mod_file = open(err_path+err_mod)
                mod = np.genfromtxt(err_mod_file)
                arg_sort = np.argsort(mod[:,1])
                val = mod[:,2][arg_sort]
                self.phase_bin = list(dict.fromkeys(mod[:,0]))
                lcf = LC_Fitter(model_name=self.model_name, sample='csp', sad_path=self.sad_path,
                                modelcov=True, qual_crit=True, version_sug=str(self.v),modeldir = self.modeldir,
                                mod_errfile = err_path+err_mod, sub_sample=val_sample, 
                                filter_drop_csp = ['cspu'], t0_fix=True, csp_file=csp_file)
                lcf.fit_sample()
                self.param_sug_path = 'param_errmod_val_'+err_mod+'.pkl'
                lcf.write_result(specific_file=self.output_path+'/'+self.param_sug_path)
                self.res_dic = {}
                try:
                    self.dic =  pkl.load(open(self.output_path+'/'+self.param_sug_path))
                except:
                    self.dic = pkl.load(open(self.output_path+'/'+self.param_sug_path,
                                             'rb'), encoding='latin1')
                self.dic = self.dic['data']
                for sn_name in self.dic.keys():
                    self.res_dic[sn_name] = {}
                    self.res_dic[sn_name]['res'], self.res_dic[sn_name]['band'], self.res_dic[sn_name]['p'] = self.residuals(sn_name)
                print(err_mod)
                likelihood = self.likelihood(val)
                if likelihood < best_likelihood:
                    best_likelihood = likelihood
                    best_node = err_mod
        print('The best error model is %s'%str(best_node))


    
                
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
        self._migrad_output_ = self.minuit.migrad(ncall=20000)
        
        if self._migrad_output_[0]["is_valid"] is False:
            print("migrad is not valid")
            
            
        self.resultsfit = np.asarray([self.minuit.values[k]
                              for k in self.freeparameters])

   
class error_mod_gp(object):
    
    def __init__(self, output_path='../../sugar_analysis_data/err_mod_training/',
                 modeldir='../../sugar_model/', sad_path = '../../',
                 bands=['cspb', 'cspg', 'cspv', 'cspr', 'cspi']):
        register_SUGAR(modeldir=modeldir, mod_errfile=modeldir+'model_err_sug.dat' ,version='0.0')
        self.bands = bands 
        self.sys = sncosmo.get_magsystem('csp')
        self.output_path = output_path
        self.modeldir = modeldir
        self.source = sncosmo.get_source('sugar', version='0.0')
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, 
                          effects=[dust], 
                          effect_names=['mw'], 
                          effect_frames=['obs']) 
        
        
    def residuals(self, sn_name):
        Filtre = self.dic[sn_name]['res']['data_mask'] 
        self.model.set(z=float(self.dic[sn_name]['res']['parameters'][0]))
        self.model.set(mwebv=self.dic[sn_name]['mwebv'])
        self.model.set(q1=float(self.dic[sn_name]['res']['parameters'][3]))
        self.model.set(q2=float(self.dic[sn_name]['res']['parameters'][4]))
        self.model.set(q3=float(self.dic[sn_name]['res']['parameters'][5]))
        self.model.set(A=float(self.dic[sn_name]['res']['parameters'][6]))
        self.model.set(Xgr=float(self.dic[sn_name]['res']['parameters'][2]))
        self.model.set(t0=float(self.dic[sn_name]['res']['parameters'][1]))
        band = self.dic[sn_name]['data_table']['band'][Filtre]
        time_obs = self.dic[sn_name]['data_table']['time'][Filtre]
        phase = (time_obs- self.dic[sn_name]['res']['parameters'][1]) / (1 + self.dic[sn_name]['res']['parameters'][0])  
        data_flux = self.dic[sn_name]['data_table']['flux'][Filtre]
        data_magerr = self.dic[sn_name]['data_table']['fluxerr'][Filtre]*1.0857362047581294/data_flux
        data_mag = np.zeros_like(data_flux)
        for j, flux in enumerate(data_flux):
            data_mag[j] = self.sys.band_flux_to_mag(flux, band[j])
        model_mag = self.model.bandmag(band, 'csp', time_obs)
        residuals = data_mag - model_mag
        return residuals, band, phase, data_magerr
    
    def build_kernel(self):
        return 0.04 * kernels.ExpSquaredKernel(6.0)

    # Define the objective function (negative log-likelihood in this case).
    def nll(self, p, y):
        self.gp.set_parameter_vector(p)
        ll = self.gp.log_likelihood(y)
        return -ll if np.isfinite(ll) else 1e25
    
    # And the gradient of the objective function.
    def grad_nll(self, p, y):
        self.gp.set_parameter_vector(p)
        return -self.gp.grad_log_likelihood(y)      
    
    def build_gp(self):
        lcf = LC_Fitter(model_name='sugar', sample='csp', sad_path=self.sad_path,
                        modelcov=False, qual_crit=True, version_sug=str(self.v),
                        modeldir = self.modeldir, sub_sample=self.training_sample,
                       filter_drop_csp = ['cspu'])
        lcf.fit_sample()
        self.param_sug_path = 'param_errmod_0it_4node.pkl'
        lcf.write_result(specific_file=self.output_path+self.param_sug_path)
        
        try:
            self.dic =  pkl.load(open(self.output_path+self.model_name+'/'+self.param_sug_path))
        except:
            self.dic = pkl.load(open(self.output_path+self.model_name+'/'+self.param_sug_path,
                                     'rb'), encoding='latin1')
        self.rids = read_input_data_SNf(res_dict_path=self.output_path+self.param_sug_path)
        self.dic = self.rids.delete_fit_fail(self.dic['data'])
        
        #Build gp object
        kernel = self.build_kernel()
        self.res_band = {}
        self.hyper = {}
        self.nb_point_per_band = {}      
        for band in self.bands:
            self.gp = george.GP(kernel, mean=0)
            self.nb_point_per_band[band] = 0.
            self.res_dic = {}
            self.res_band[band] = {}
            
            self.res_band[band]['res'] = []
            self.res_band[band]['phase'] = []
            self.res_band[band]['err'] = []
            for sn_name in self.dic.keys():
                self.res_dic[sn_name] = {}
                self.res_dic[sn_name]['res'], self.res_dic[sn_name]['band'], self.res_dic[sn_name]['phase'], self.res_dic[sn_name]['err'] = self.residuals(sn_name)
                
                Filtre = np.ma.masked_array(self.res_dic[sn_name]['band'], mask=(self.res_dic[sn_name]['band'] == band))
                self.res_band[band]['res'] += list(self.res_dic[sn_name]['res'][Filtre.mask])
                self.res_band[band]['phase'] += list(self.res_dic[sn_name]['phase'][Filtre.mask])
                self.res_band[band]['err'] += list(self.res_dic[sn_name]['err'][Filtre.mask])
                #self.nb_point_per_band[bands] += len(self.res_dic[sn_name]) 
            
            self.gp.compute(np.array(self.res_band[band]['phase']), np.array(self.res_band[band]['err']))
            # minimize to find the best-fit parameters
            initial_params = self.gp.get_parameter_vector()
            bounds = self.gp.get_parameter_bounds()
            try :
                self.hyper[band] = minimize(self.nll, initial_params, jac=self.grad_nll, 
                             method="L-BFGS-B", bounds=bounds, args=(np.array(self.res_band[band]['res'])))
                print(self.hyper[band])
            except:
                print(band+' failed')