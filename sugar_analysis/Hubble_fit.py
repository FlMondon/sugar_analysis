#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 15:56:16 2018

@author: florian
"""

import sncosmo
import builtins
import numpy as np
import cosmo_tools as ct
from scipy.linalg import block_diag
from numpy.linalg import inv
import cPickle as pkl
import iminuit as minuit
from scipy import optimize
import math_toolbox as math
import copy

output_path = '../../'


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

class read_input_data_SNf(object):
    """
    """
    def __init__(self, res_dict_path=output_path, model_name='sugar', select=None, standard_SNf='mb'):
        self.model_name = model_name
        dic =  pkl.load(open(res_dict_path+'sugar_analysis_data/resfitlc_SNf_'+self.model_name+'.pkl'))
        self.dic_res = dic['data']
        self.t0_fix = dic['info']['t0 fix']
        self.select = select
        self.standard_SNf = standard_SNf
    
    def delete_fit_fail(self):
        """
        """
        self.fit_fail = []
        dic_del = {}
        for sn_name in self.dic_res.keys():
            if self.dic_res[sn_name]['res'] != 'fit fail':
                dic_del[sn_name] = self.dic_res[sn_name]
            else:
                self.fit_fail.append(sn_name)
        self.dic_res = dic_del 
    
    def selection(self):
        
        if self.select is not None:
            dic = {}
            sn_out_selec = []
            for sn_name in self.select:
                if sn_name in self.dic_res.keys():
                    dic[sn_name] = self.dic_res[sn_name]
                else:
                    sn_out_selec.append(sn_name)
            self.sn_out_selec = sn_out_selec
            self.dic_res = dic
                   
    def trans_data_form(self):
        self.delete_fit_fail()
        self.selection()
        list_param = []
        list_cov = []
        list_zhl = []
        list_zcmb = []
        list_zerr = []
        for sn_name in self.dic_res.keys():
            
            list_param.append(list(self.dic_res[sn_name]['res']['parameters'][2:-2]))
            if self.t0_fix:
                list_cov.append(self.dic_res[sn_name]['res']['covariance'])
            else:
                list_cov.append(np.delete(np.delete(self.dic_res[sn_name]['res']['covariance'], 0, 0),0,1))
            list_zhl.append(self.dic_res[sn_name]['zhl'])
            list_zcmb.append(self.dic_res[sn_name]['zcmb'])
            list_zerr.append(self.dic_res[sn_name]['zerr'])
        self.zhl = np.array(list_zhl)
        self.zcmb = np.array(list_zcmb)
        self.zerr = np.array(list_zerr)
        self.params = np.array(list_param)
        self.cov = np.array(list_cov)
        
    def Xgr_to_Mgr(self):
        for j in range(len(self.cov[:,0])):
            self.cov[j,0,:] = self.cov[j,0,:]* (1./(-2.5*np.log(10)*self.params[j,0]))
            self.cov[j,:,0] = self.cov[j,:,0]* (1./(-2.5*np.log(10)*self.params[j,0]))
        self.params[:,0] = -2.5*np.log10(self.params[:,0])

    def  Xgr_to_mb(self):
        builtins.builtins_jla_bandpasses()
        builtins.mag_sys_jla()
        # Calculation of m_b
        builtins.register_SUGAR()
        source = sncosmo.get_source('sugar')
        model = sncosmo.Model(source=source)
        self.mb = np.zeros(len(self.dic_res.keys()))
        cov_copy = copy.deepcopy(self.cov)
        for i, sn_name in enumerate (self.dic_res.keys()):
            parameters = self.dic_res[sn_name]['res']['parameters']
            model.set(t0=parameters[1], Xgr=parameters[2], q1=parameters[3],
                       q2=parameters[4], q3=parameters[5], A=parameters[6])
            self.mb[i] = model.bandmag('jla_STANDARD::B', 'jla_VEGA2_mb', [parameters[1]]) 
            h = 10**-9
            #dmbxgr
            model.set(t0=parameters[1], Xgr=parameters[2]+h, q1=parameters[3],
                       q2=parameters[4], q3=parameters[5], A=parameters[6])
            dmbxgr = (model.bandmag('jla_STANDARD::B', 'jla_VEGA2_mb', [parameters[1]]) - self.mb[i])/h
            #dmbq1
            model.set(t0=parameters[1], Xgr=parameters[2], q1=parameters[3]+h,
                       q2=parameters[4], q3=parameters[5], A=parameters[6])
            dmbq1 = (model.bandmag('jla_STANDARD::B', 'jla_VEGA2_mb', [parameters[1]]) - self.mb[i])/h
            #dmbq2
            model.set(t0=parameters[1], Xgr=parameters[2], q1=parameters[3],
                       q2=parameters[4]+h, q3=parameters[5], A=parameters[6])
            dmbq2 = (model.bandmag('jla_STANDARD::B', 'jla_VEGA2_mb', [parameters[1]]) - self.mb[i])/h
            #dmbq3
            model.set(t0=parameters[1], Xgr=parameters[2], q1=parameters[3],
                       q2=parameters[4], q3=parameters[5]+h, A=parameters[6])
            dmbq3 = (model.bandmag('jla_STANDARD::B', 'jla_VEGA2_mb', [parameters[1]]) - self.mb[i])/h
            #dmbA
            model.set(t0=parameters[1], Xgr=parameters[2], q1=parameters[3],
                       q2=parameters[4], q3=parameters[5], A=parameters[6]+h)            
            dmbA = (model.bandmag('jla_STANDARD::B', 'jla_VEGA2_mb', [parameters[1]]) - self.mb[i])/h
            
            dray = np.array([dmbxgr, dmbq1, dmbq2, dmbq3, dmbA])
            for j in range(len(parameters[2:-2])):
                if j ==0:
                    self.cov[i,0,j] = np.dot(np.dot(dray.T,cov_copy[i]),dray)
                else:
                    self.cov[i,0,j] = np.dot(cov_copy[i,j,:],dray)
                    self.cov[i,j,0] = self.cov[i,0,j]
        self.params[:,0] = self.mb
    def build_HD_data(self):
        self.trans_data_form()
        if self.standard_SNf == 'Mgr':
            self.Xgr_to_Mgr()
        elif self.standard_SNf == 'mb':
            self.Xgr_to_mb()
        else:
            raise ValueError('standard_SNf have to be Mgr or mb')
        self.cov = list(self.cov)
        self.cov = block_diag(*self.cov)
    
          

def get_hubblefit(x, cov_x, zhl, zcmb, zerr, PARAM_NAME=np.asarray(['alpha1','alpha2',"alpha3","beta","delta"])):
    """
    Parameters
    ----------
    x: type
        infor
    
    cov_x ....
    
    
    parameters: list of int to specify if you want to remove some parameters in the fit in PARAM_NAME. 
                By default is None and this take all parameters. If you don't want to put a correction use 
                parameters=[]
        example: 
            for usual salt x1 c correction if you want only color correction, use parameters = [2]. 
            Warning if you have more than 5 parameters, add parameters in PARAM_NAME
    """
    n_corr =  np.shape(x)[1]-1
    
    class hubble_fit_case(Hubble_fit):
        freeparameters = ["Mb"]+PARAM_NAME[:n_corr].tolist()

    h = hubble_fit_case(x, cov_x, zhl, zcmb, zerr)
    return h


class Hubble_fit(object):
    """
    """
    
    def __new__(cls,*arg,**kwargs):
        """ Upgrade of the New function to enable the
        the _minuit_ black magic
        """
        obj = super(Hubble_fit,cls).__new__(cls)
        
        exec "@make_method(Hubble_fit)\n"+\
             "def _minuit_chi2_(self,%s): \n"%(", ".join(obj.freeparameters))+\
             "    parameters = %s \n"%(", ".join(obj.freeparameters))+\
             "    return self.get_chi2(parameters)\n"


        return obj
    
    def __init__(self, X, cov_X, zhl, zcmb, zerr, guess=None):
        
        self.variable = X
        self.cov = cov_X
        self.zcmb = zcmb
        self.zhl = zhl
        self.zerr = zerr
        self.dmz = 5/np.log(10) * np.sqrt(self.zerr**2 + 0.001**2) / self.zcmb
        self.dof = len(X)-len(self.freeparameters)  
        
    

    

    def distance_modulus(self, params):
        """
        (mb + alpha * v1 + beta * v2 .....) - Mb
        """
        return  np.sum(np.concatenate([[1],params[1:]]).T * self.variable, axis=1) - params[0]
                
    def get_chi2(self, params):
        """
        """
        self.Cmu = np.zeros_like(self.cov[::len(params),::len(params)])
        pcorr = np.concatenate([[1],params[1:]])
        for i, coef1 in enumerate(pcorr):
            for j, coef2 in enumerate(pcorr):
                self.Cmu += (coef1 * coef2) * self.cov[i::len(params),j::len(params)]  
                
                
        self.Cmu[np.diag_indices_from(self.Cmu)] += self.sig_int**2 + self.dmz**2 
        self.C = inv(self.Cmu)
        L = self.distance_modulus(params) - ct.distance_modulus_th(self.zcmb, self.zhl)
        self.residuals = L
        self.var = np.diag(self.Cmu)
        return np.dot(L, np.dot(self.C,L))        


    def fit_intrinsic(self, intrinsic_guess=0.1):
        """ Get the most optimal intrinsic dispersion given the current fitted standardization parameters. 
        
        The optimal intrinsic magnitude dispersion is the value that has to be added in quadrature to 
        the magnitude errors such that the chi2/dof is 1.
        Returns
        -------
        float (intrinsic magnitude dispersion)
        """
        def get_intrinsic_chi2dof(intrinsic):
            self.sig_int = intrinsic
            return np.abs(self.get_chi2(self.resultsfit) / self.dof -1)
        
        return optimize.fmin(get_intrinsic_chi2dof,
                             intrinsic_guess, disp=0)[0]    

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
        def _test_it_(k,info):
            param = k.split(info)[0]
            if param not in self.freeparameters:
                raise ValueError("Unknown parameter %s"%param)

        self.param_input = {}
        # -- what you hard coded
        for name in self.freeparameters:
            for info in ["_guess","_fixed","_boundaries"]:
                if hasattr(self, name+info):
                    self.param_input[name+info] = eval("self.%s"%(name+info))
                    
        # -- Then, whatever you gave
        for k,v in kwargs.items():
            if "_guess" in k:
                _test_it_(k,"_guess")
            elif "_fixed" in k:
                _test_it_(k,"_fixed")
            elif "_boundaries" in k:
                _test_it_(k,"_boundaries")
            else:
                raise ValueError("I am not able to parse %s ; not _guess, _fixed nor _boundaries"%k)
            self.param_input[k] = v

        # -- Finally if no values have been set, let's do it
        for name in self.freeparameters:
            if name+"_guess" not in self.param_input.keys():
                self.param_input[name+"_guess"] = 0.
            if name+"_fixed" not in self.param_input.keys():
                self.param_input[name+"_fixed"] = False
            if name+"_boundaries" not in self.param_input.keys():
                self.param_input[name+"_boundaries"] = [None,None]
    
    def fit(self, fit_intrinsic=True, **kwargs):
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
        self._loopcount = 0
        self.sig_int = 0.
        self.setup_guesses(**kwargs)
        
        self.first_iter = self._fit_minuit_()
        # - Intrinsic disposerion Fit?
        if fit_intrinsic:
            while (np.abs(self.chi2_per_dof - 1) > 0.001 and self._loopcount < 30):
                self.sig_int =  self.fit_intrinsic(np.sqrt(np.mean(self.var))*2. / self.chi2_per_dof)
                self._fit_minuit_()
                self._loopcount += 1
                
        # - Final steps      
        return self._fit_readout_()

        
    def _setup_minuit_(self):
        """
        """
        # == Minuit Keys == #
        minuit_kwargs = {}
        for param in self.freeparameters:
            minuit_kwargs[param] = self.param_input["%s_guess"%param]


        self.minuit = minuit.Minuit(self._minuit_chi2_, **minuit_kwargs)
    
    def _fit_minuit_(self):
        """
        """
        self._setup_minuit_()
        self._migrad_output_ = self.minuit.migrad()
        
        if self._migrad_output_[0]["is_valid"] is False:
            print("migrad is not valid")
            
            
        self.resultsfit = np.asarray([self.minuit.values[k]
                              for k in self.freeparameters])
        self.chi2_per_dof = self.minuit.fval/self.dof
        
    def _fit_readout_(self):
        """ Computes the main numbers """
        return math.comp_rms(self.residuals, self.dof, err=True, variance=self.var) , self.sig_int  
    

            
