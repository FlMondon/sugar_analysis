#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 10:40:00 2019

@author: florian
"""

import sncosmo
import numpy as np
import os
from iminuit import Minuit
import pickle as pkl
from .constant import t_min_sug, t_max_sug 
from .fitting_lc import LC_Fitter
from .Hubble_fit import read_input_data_SNf 
from .load_sugext import register_SUGAR_ext

class train_transext(object):
    
    def __init__(self, modelcov=True, sad_path = '../../',
                 modeldir='../../sugar_model/',
                 res_path_salt2='../../sugar_analysis_data/resfitlc_snls_salt2.pkl',
                 mod_errfile='../../sugar_analysis_data/err_mod_training/sugar/err_mod_4node.dat',
                 sn_subsample='SNLS'): 
            
            
            self.res_path_salt2 = res_path_salt2
            self.sad_path = sad_path
            self.modeldir = modeldir
            self.mod_errfile = mod_errfile
            self.modelcov = modelcov
            self.sn_subsample = sn_subsample
            register_SUGAR_ext(modeldir=self.modeldir, mod_errfile=self.mod_errfile,
                       transmat_path='trans_matrix_init.pkl', version='0.0')
            self.lc_fitter = LC_Fitter(model_name='sugar_ext', sample='jla', param_sug =['q1', 'q3'],
                                       sub_sample=self.sn_subsample,
                                       modelcov=self.modelcov, version_sug='0.0')
    
    def print_mean_chi(self):
        chi_salt2 = []
        chi_sugar = []
        for sn_name in self.rids_salt2.dic_res.keys():
            chi_salt2.append(self.rids_salt2.dic_res[sn_name]['res']['chisq']/self.rids_salt2.dic_res[sn_name]['res']['ndof'])
            chi_sugar.append(self.rids_sugar.dic_res[sn_name]['res']['chisq']/self.rids_sugar.dic_res[sn_name]['res']['ndof'])
        print('mean chi2 per dof salt2: ', np.mean(chi_salt2), 'mean chi2 per dof sugar it n-1: ', np.mean(chi_sugar))   
        
    def train(self):
            self.lc_fitter.fit_sample()
            res_path = os.path.join(self.sad_path,'sugar_analysis_data/extension_training/resfitlc_snls_ext_0.pkl')
            self.lc_fitter.write_result(specific_file= res_path)
            n = 1
            while n <= 3:
                pksalt2 = pkl.load(open(self.res_path_salt2,'rb'))
                salt2_sn_list = [sn for sn in pksalt2['data'] if pksalt2['data'][sn]['res']!='fit fail']
                pksug = pkl.load(open(res_path,'rb'))         
                sug_sn_list = [sn for sn in pksug['data'] if pksug['data'][sn]['res']!='fit fail']
                common_sn_list = list(set(sug_sn_list) & set(salt2_sn_list))
                self.rids_sugar = read_input_data_SNf(select=common_sn_list, res_dict_path=res_path, standard='Mgr')
                self.rids_sugar.build_HD_data()
                self.rids_salt2 = read_input_data_SNf(select=self.rids_sugar.sn_name, res_dict_path=self.res_path_salt2, model_name='salt2', standard='log10_x0')
                self.rids_salt2.build_HD_data()
                sug_params = self.rids_sugar.params
                salt2_params = self.rids_salt2.params
                nlines,ncol = np.shape(sug_params)
                Q = np.ones((nlines, ncol+1))
                Q[:,1:]=sug_params
                sQQT = Q.T.dot(Q)
                sQXt = Q.T.dot(salt2_params)
                A = np.linalg.inv(sQQT).dot(sQXt)
                pkl.dump(A, open('../../sugar_model/trans_matrix'+str(n)+'.pkl', 'wb'))
                register_SUGAR_ext(modeldir=self.modeldir, mod_errfile=self.mod_errfile,
                           transmat_path='trans_matrix'+str(n)+'.pkl', version=str(n)+'.0')
                self.lc_fitter = LC_Fitter(model_name='sugar_ext', sample='jla', param_sug =['q1', 'q3'],
                                           modelcov=self.modelcov, sub_sample=self.sn_subsample, version_sug=str(n)+'.0')       
                self.lc_fitter.fit_sample()
                res_path = os.path.join(self.sad_path,'sugar_analysis_data/extension_training/resfitlc_snls_ext_'+str(n)+'.pkl')
                self.lc_fitter.write_result(specific_file= res_path)
                n += 1
                self.print_mean_chi()

    