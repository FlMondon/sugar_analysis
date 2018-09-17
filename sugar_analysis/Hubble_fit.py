#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 12 15:56:16 2018

@author: florian
"""

import numpy as np
import constant as cst
from scipy import integrate
from scipy.linalg import block_diag
from numpy.linalg import inv
import cPickle as pkl
output_path = '../../'



class read_input_data_SNf(object):
    """
    """
    def __init__(self, res_dict_path=output_path, model_name='sugar'):
        self.model_name = model_name
        self.dic_res =  pkl.load(open(res_dict_path+'sugar_analysis_data/resfitlc_SNf_'+self.model_name+'.pkl'))
    
    def delete_fit_fail(self):
        """
        """
        self.fit_fail = []
        self.dic_del = {}
        for sn_name in self.dic_res.keys():
            if self.dic_res[sn_name]['res'] != 'fit fail':
                self.dic_del[sn_name] = self.dic_res[sn_name]['res']
                   
    def trans_data_form(self):
        self.delete_fit_fail()
        list_param = []
        list_cov = []
        list_zhl = []
        list_zcmb = []
        list_zerr = []
        for sn_name in self.dic_del.keys():
            
            list_param.append(list(self.dic_del[sn_name]['parameters'][2:-2]))
            list_cov.append(np.delete(np.delete(self.dic_del[sn_name]['covariance'], 0, 0),0,1))
            list_zhl.append(self.dic_del[sn_name]['zhl'])
            list_zcmb.append(self.dic_del[sn_name]['zcmb'])
            list_zerr.append(self.dic_del[sn_name]['zerr'])
        self.zhl = np.array(list_zhl)
        self.zcmb = np.array(list_zcmb)
        self.zerr = np.array(list_zerr)
        self.params = np.array(list_param)
        self.cov = np.array(list_cov)
        
    def Xgr_to_Mgr(self):
        self.trans_data_form()
        for j in range(len(self.cov[:,0])):
            self.cov[j,0,:] = self.cov[j,0,:]* (1./(2.5*np.log(10)*self.params[j,0]))
            self.cov[j,:,0] = self.cov[j,:,0]* (1./(2.5*np.log(10)*self.params[j,0]))
        self.params[:,0] = -2.5*np.log10(self.params[:,0])

    def build_HD_data(self):
        self.Xgr_to_Mgr()
        self.cov = list(self.cov)
        self.cov = block_diag(*self.cov)
    
class Hubble_fit_sugar(object):
    """
    """
    def __init__(self, X, cov_X, zhl, zcmb, zerr, guess=None):
        
        self.params = X
        self.cov = cov_X
        self.zhl = zhl
        self.zcmb = zcmb
        self.zerr = zerr
        self.dmz = 5/np.log(10) * np.sqrt(self.zerr**2 + 0.001**2) / self.zcmb
        self.dof = np.shape(X)[0] - np.shape(X)[1]
        self.nb_params = np.shape(X)[1]
        
        if guess==None:
            self.guess = np.zeros(self.nb_params)
        else:
            self.guess = guess

    # ------------------ #
    #   Cosmology        #
    # ------------------ #

    def int_cosmo(self, z, Omega_M=0.3):   
        """
        """
        return 1./np.sqrt(Omega_M*(1+z)**3+(1.-Omega_M))
        
    def luminosity_distance(self):
        """
        """
        if type(self.zcmb)==np.ndarray:
            integr = np.zeros_like(self.zcmb)
            for i in range(len(self.zcmb)):
                integr[i]=integrate.quad(self.int_cosmo, 0, self.zcmb[i])[0]
        else:
            integr = integrate.quad(self.int_cosmo, 0, self.zcmb)[0]
    
        return (1+self.zhl)*(cst.CLIGHT/cst.H0)*integr
 
    def distance_modulus_th(self):
        """
        """
        return 5.*np.log(self.luminosity_distance())/np.log(10.)-5.
    

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
        L = self.distance_modulus(params) - self.distance_modulus_th()
        self.residuals = L
        self.var = np.diag(self.Cmu)
        return np.dot(L, np.dot(self.C,L))        