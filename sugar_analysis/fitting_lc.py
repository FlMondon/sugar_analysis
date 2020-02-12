#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 10:33:31 2018

@author: florian
"""

import sncosmo
import copy
import numpy as np
import os
import pickle as pkl
   
from .builtins import register_SNf_bands_width, mag_sys_SNF_width,  builtins_jla_bandpasses, mag_sys_jla
#from .load_sugar import register_SUGAR
from .read_jla import read_lc_jla
from .cosmo_tools import distance_modulus_th
from .constant import t_min_sug, t_max_sug, t_min_salt2, t_max_salt2
from .data_table import build_data
from .data_table import read_csp
from .load_salt2_newerr import register_salt2_newerr
from .load_sugext import register_SUGAR_ext
   
class LC_Fitter(object):
    
    def __init__(self, model_name='sugar', sample='SNf', data=None, 
                 t0_fix=False, sub_sample=None, modelcov=True,
                 filters=['new_fU_10','fB_10','fV_10','fR_10','new_fI_10'], 
                 filter_drop_csp = None, sad_path = '../../',
                 modeldir='../../sugar_model/', mod_errfile='model_err_sug.dat',
                 width=10, param_sug =['q1', 'q2', 'q3'], qual_crit=False, 
                 version_sug='1.0'):
        
        
        self.fitting_sample =  False
        self.model_name = model_name
        self.t0_fix = t0_fix
        self.sad_path = sad_path
        self.version_sug = version_sug
        self.filter_drop_csp = filter_drop_csp
        self.dic_res = None
        self.sub_sample = sub_sample
        self.sample = sample
        self.qual_crit = qual_crit
        self.modeldir = modeldir
        self.modelcov = modelcov
        self.filters = filters
        self.strfilters = str(len(filters))+filters[0]+filters[len(filters)-1]
#        register_SUGAR(modeldir=self.modeldir, mod_errfile=mod_errfile,
#                       version=self.version_sug)
        if model_name == 'salt2_newerr':
            register_salt2_newerr(modeldir=modeldir, mod_errfile=mod_errfile,
                       version=self.version_sug)
        self.param_sug = param_sug

                
        
        if self.sample=='SNf':
            self.bd = build_data(sad_path=self.sad_path, modeldir=self.modeldir)
            self.meta = pkl.load(open(self.sad_path+'sugar_analysis_data/SNF-0203-CABALLO2/META.pkl','rb'), encoding="latin1")
            self.data = []
            if sub_sample==None:
                for sn_name in self.meta.keys():
                    if self.meta[sn_name]['idr.subset'] == 'training' or self.meta[sn_name]['idr.subset'] == 'validation':
                        self.data.append(sn_name)
            else :
                for sn_name in self.sub_sample:
                    self.data.append(sn_name)

            self.errorscale = True
            self.width= width
            register_SNf_bands_width(width=self.width, sad_path=self.sad_path)
            mag_sys_SNF_width(width=self.width, sad_path=self.sad_path)
            
        elif sample=='jla':
            builtins_jla_bandpasses(sad_path=self.sad_path)
            mag_sys_jla(sad_path=self.sad_path)
            self.dic_jla_zbias = {}
            jla_file = np.loadtxt(self.sad_path+'sugar_analysis_data/data/jla_data/jla_lcparams.txt',dtype='str')
            datos = os.listdir(self.sad_path+'sugar_analysis_data/data/jla_data/jla_light_curves/')
            self.data = []
            for sn_name in datos:
                if sub_sample==None:
                    if sn_name.startswith('lc-'):
                        self.data.append(sn_name)
                elif sub_sample=='SDSS':
                    if sn_name.startswith('lc-SDSS'):
                        self.data.append(sn_name)      
                elif sub_sample=='SNLS':
                    if sn_name.startswith('lc-0'):
                        self.data.append(sn_name)                     
                elif sub_sample=='Nearby':
                    if sn_name.startswith('lc-sn'):
                        self.data.append(sn_name) 
                elif sub_sample=='HST':
                    HST_name = ['lc-Vilas.list','lc-Torngasek.list',
                                'lc-Ombo.list', 'lc-Patuxent.list',
                                'lc-Lancaster.list', 'lc-Gabi.list',
                                'lc-Gabi.list','lc-Eagle.list',
                                'lc-Aphrodite.list']
                    if sn_name in HST_name:
                        self.data.append(sn_name)              
                else:
                    if sn_name in self.sub_sample:
                        self.data.append(sn_name)   
            for line in jla_file:
                if 'lc-' + line[0] + '.list' in self.data:
                    self.dic_jla_zbias['lc-' + line[0] + '.list'] = float(line[1]), float(line[3]), float(line[20])
        
        elif sample=='csp':
            builtins_jla_bandpasses(sad_path=self.sad_path)
            mag_sys_jla(sad_path=self.sad_path)
            self.rcsp = read_csp(sad_path=self.sad_path)
            datos = os.listdir(self.sad_path+'sugar_analysis_data/DR3/')
            self.data = []
            csp_type_file = open(self.sad_path+'sugar_analysis_data/DR3/krisciunas17_table2.org')
            csp_type_lines = csp_type_file.readlines()
            dic_csp_type={}
            for line in csp_type_lines:
                line = line.split()
                dic_csp_type['SN'+line[1] +'_snpy.txt'] =line[5]
                
            for sn_name in datos:
                if sub_sample==None:
                    if sn_name.startswith('SN2') and dic_csp_type[sn_name]=='normal':
                        self.data.append(sn_name)  
                else:
                    if sn_name in self.sub_sample and dic_csp_type[sn_name]=='normal':
                        self.data.append(sn_name) 
            
            self.dic_csp_radec = {}
            csp_file = open(self.sad_path+'sugar_analysis_data/DR3/tab1.dat')
            csp_line = csp_file.readlines()
            for line in csp_line:
                line = line.split()
                if 'SN'+line[0] + '_snpy.txt' in self.data:
                    self.dic_csp_radec['SN'+line[0] +'_snpy.txt'] =line[1], line[2]
            

            
        else:
            self.width= width
            register_SNf_bands_width(width=self.width)
            mag_sys_SNF_width(width=self.width)
            Warning('Mag sys and band used in your data have to be register in sncosmo our in builtins')                
            if type(data) is not dict and data is not None:
                print (type(data))
                raise ValueError('data have to be a dictionary or None')
            elif data is not None:
                self.data = pkl.load(open(data))
            else:
                Warning('Data is not defined you can only fit a sigle SNIa with an astropy Table')
    def sn_data(self, sn_name):
       
        if self.sample == 'SNf':
            self.zhl, self.zcmb, self.zerr, self.mwebv, self.daymax = self.meta[sn_name]['host.zhelio'], self.meta[sn_name]['host.zcmb'], self.meta[sn_name]['host.zhelio.err'], self.meta[sn_name]['target.mwebv'], self.meta[sn_name]['salt2.DayMax']
            self.table_sn = self.bd.build_Astropy_Table(sn_name, band_used=self.filters)
        
        elif self.sample == 'jla':
            self.head, self.table_sn = read_lc_jla(sn_name, sad_path=self.sad_path, model=self.model_name)
            self.zhl, self.zcmb, self.zerr, self.biascor, self.mwebv = self.head['@Z_HELIO'], self.dic_jla_zbias[sn_name][0], self.dic_jla_zbias[sn_name][1], self.dic_jla_zbias[sn_name][2], self.head['@MWEBV']

        elif self.sample == 'csp':
            self.table_sn, self.zhl = self.rcsp.build_csp_table(sn_name, drop=self.filter_drop_csp)
            self.mwebv, self.zcmb = self.rcsp.get_mwebv(self.dic_csp_radec[sn_name][0], self.dic_csp_radec[sn_name][1], self.zhl)
            self.zerr = 0.
        else:
            self.table_sn = self.data[sn_name]['table']
            self.zhl = self.data[sn_name]['zhl']
            self.zerr = self.data[sn_name]['zerr']
            self.zcmb = self.data[sn_name]['zcmb']
            self.mwebv = self.data[sn_name]['mwebv']
            if self.t0_fix== True:
                try : 
                    self.daymax = self.data[sn_name]['daymax']
                except:
                    raise ValueError('t0_fix is True: need daymax in the dictionary data')
            
    
    def fit_lc_sugar(self, data, zhl=None, zcmb=None, mwebv=None):
        """
        """
        if not self.fitting_sample :
            self.source = sncosmo.get_source('sugar', version=self.version_sug)
            dust = sncosmo.CCM89Dust()
            self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs']) 
            self.mwebv = mwebv
            self.zhl = zhl
            #If you don't have zcmb we used only zhl to compute the distance
            #modulus for the guess of Xgr
            if zcmb == None:
                self.zcmb = zhl
            else:
                self.zcmb = zcmb
            if mwebv == None:
                self.model.set(mwebv=0.)
            else: 
                self.model.set(mwebv=self.mwebv)
            if zhl==None:
                self.model.set(z=0.)
                Xgr_init = 1.e-15
            else:
                self.model.set(z=self.zhl)
                Xgr_init = 10**(-0.4*(distance_modulus_th(self.zcmb,self.zhl)))
        else:
            self.model.set(z=self.zhl)
            self.model.set(mwebv=self.mwebv)
            if self.zcmb==None:
                Xgr_init = 10**(-0.4*(distance_modulus_th(self.zhl,self.zhl)))
            else:
                Xgr_init = 10**(-0.4*(distance_modulus_th(self.zcmb,self.zhl)))
        self.model.set(Xgr=Xgr_init)
        if self.t0_fix:
            self.model.set(t0=self.daymax)
            res, fitted_model = sncosmo.fit_lc(data, 
                                               self.model, 
                                               self.param_sug+['A', 'Xgr'], 
                                               modelcov  = self.modelcov)
        else:
            res, fitted_model = sncosmo.fit_lc(data, self.model, 
                                            ['t0']+self.param_sug+['A', 'Xgr'], 
                                            modelcov  = self.modelcov,
                                               phase_range=(t_min_sug, 
                                                            t_max_sug))
            
        return res, fitted_model 
    
    

        
    def fit_lc_salt2(self, data, zhl=None,  mwebv=None ):
        """
        """
        if not self.fitting_sample :
            if self.model_name == 'salt2' :
                self.source = sncosmo.get_source(self.model_name, version='2.4')
            else:
                self.source = sncosmo.get_source(self.model_name, version=self.version_sug)
            dust = sncosmo.CCM89Dust()
            self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs'])     
            if zhl==None:
                try:
                    self.model.set(z=self.zhl)
                except:
                    self.model.set(z=0.)
            else:
                self.model.set(z=zhl)
            if mwebv==None:
                self.model.set(mwebv=0.)
            else:
                self.model.set(mwebv=mwebv)
        else:
            self.model.set(z=self.zhl)
            self.model.set(mwebv=self.mwebv)
        if self.t0_fix:
            self.model.set(t0=self.daymax)
            res, fitted_model = sncosmo.fit_lc(data, self.model, 
                                                   ['x1',
                                                    'c', 
                                                    'x0'], 
                                                   modelcov  = self.modelcov)
        else:
            res, fitted_model = sncosmo.fit_lc(data, self.model, 
                                                   ['t0',
                                                    'x1',
                                                    'c', 
                                                    'x0'], 
                                                   modelcov  = self.modelcov,
                                                   phase_range=(t_min_salt2, 
                                                                t_max_salt2))
        
        return res, fitted_model
        
    def selection_new_criteria(self, data, zhl, mwebv):
        
        if len(data) >= 4.:
            if self.model_name == 'sugar' or self.model_name == 'sugar_ext':
                paran_sug_init = copy.deepcopy(self.param_sug)
                self.param_sug = []
                try:
                    zcmb = self.zcmb
                except:
                    zcmb = None
                try:
                    res, fitted_model = self.fit_lc_sugar(data, zhl=zhl, 
                                                      mwebv=mwebv,
                                                      zcmb=zcmb)
                    fit_sucess = True
                    
                except:
                    bool_sectect = False
                    fit_sucess = False
                if fit_sucess:
                    t_peak = fitted_model.parameters[1]
                    t1 = t_peak -12*(1 + self.model.get('z'))
                    t2 = t_peak + 5*(1 + self.model.get('z'))
                    t3 = t_peak + 20*(1 + self.model.get('z'))     
                    t4 = t_peak -8*(1 + self.model.get('z'))
                    t5 = t_peak + 10*(1 + self.model.get('z'))
                    
                    crit1 = False
                    crit2 = False
                    crit3 = False
                    nb_filter1 = 0
                    prev_filter = []
                    for j, time in enumerate(data['time']):
                        if time >= t1 and time <=t2:
                                crit1 = True
                        if time >= t2 and time <=t3:
                            crit2 = True
                        if time >= t4 and time <=t5:
                            if nb_filter1 >= 1 and data['band'][j] not in prev_filter:
                                crit3 = True
                            else:
                                nb_filter1 += 1 
                                prev_filter.append(data['band'][j])
                        
                    
                    bool_sectect = crit1 and crit2 and crit3
                self.param_sug = paran_sug_init
            elif self.model_name == 'salt2' or self.model_name == 'salt2_newerr' :
                if not self.fitting_sample :
                    if self.model_name == 'salt2' :
                        self.source = sncosmo.get_source(self.model_name, version='2.4')
                    else:
                        self.source = sncosmo.get_source(self.model_name, version=self.version_sug)
                    dust = sncosmo.CCM89Dust()
                    self.model = sncosmo.Model(source=self.source, 
                                      effects=[dust], 
                                      effect_names=['mw'], 
                                      effect_frames=['obs'])         
                self.model.set(mwebv=mwebv)
                self.model.set(z=zhl)
                try:
                    res, fitted_model = sncosmo.fit_lc(data, self.model, 
                                                       ['t0', 
                                                        'c', 
                                                        'x0'], 
                                                       modelcov  = True,
                                                       phase_range=(t_min_salt2, 
                                                                    t_max_salt2))
                    fit_sucess = True
                    
                except:
                    bool_sectect = False
                    fit_sucess = False
                if fit_sucess:
                    t_peak = fitted_model.parameters[1]
                    t1 = t_peak -10*(1 + self.model.get('z'))
                    t2 = t_peak + -3*(1 + self.model.get('z'))
                    t3 = t_peak + 20*(1 + self.model.get('z'))     
                    t4 = t_peak -8*(1 + self.model.get('z'))
                    t5 = t_peak + 10*(1 + self.model.get('z'))
                    
                    crit1 = False
                    crit2 = False
                    crit3 = False
                    nb_filter = 0
                    prev_filter = []
                    for j, time in enumerate(data['time']):
                        if time >= t1 and time <=t2:
                            crit1 = True
                        if time >= t2 and time <=t3:
                            crit2 = True
                        if time >= t4 and time <=t5:
                            if nb_filter >= 1 and data['band'][j] not in prev_filter:
                                crit3 = True
                            else:
                                nb_filter += 1 
                                prev_filter.append(data['band'][j])
                        
                    
                    bool_sectect = crit1 and crit2 and crit3
            else:
                raise ValueError ('model_name have to be sugar or salt2')
        else:
            bool_sectect = False
            
        return bool_sectect
    
    def fit_sample(self):
        """
        """
        self.dic_res = {}
        self.dic_res['data'] = {}
        self.dic_res['info'] = {}
        self.dic_res['info']['t0 fix'] = self.t0_fix
        self.dic_res['info']['model used'] = self.model_name
        self.dic_res['info']['sample'] = self.sample
        if self.sub_sample is not None:
            self.dic_res['info']['sub sample'] = self.sub_sample
        self.fit_fail = []
        self.fitted_model = []
        self.fitting_sample = True
        if self.model_name == 'salt2':
            self.source = sncosmo.get_source(self.model_name, version='2.4')
        else:
            self.source = sncosmo.get_source(self.model_name, version=self.version_sug)
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs'])
        for sn_name in self.data:
            self.sn_data(sn_name)
            print (sn_name)
            
            self.dic_res['data'][sn_name] = {}
#            if self.qual_crit:
#                bool_selct = self.selection_new_criteria(self.table_sn,
#                                                         self.zhl, self.mwebv)
#                if bool_selct:
#                    self.dic_res['data'][sn_name]['selection'] = True
#                else:
#                    res_sn, fitted_model_sn = 'fit fail', np.nan
#                    self.dic_res['data'][sn_name]['selection'] = False
#                    self.dic_res['data'][sn_name]['data_table'] = np.array(self.table_sn)
#                    self.dic_res['data'][sn_name]['res'] = res_sn
#                    self.dic_res['data'][sn_name]['zhl'] = self.zhl
#                    self.dic_res['data'][sn_name]['mwebv'] = self.mwebv 
#                    self.dic_res['data'][sn_name]['zerr'] = self.zerr
#                    self.dic_res['data'][sn_name]['zcmb'] = self.zcmb   
#                    continue
            if self.model_name == 'sugar' or self.model_name == 'sugar_ext':
                try:
                    res_sn, fitted_model_sn = self.fit_lc_sugar(self.table_sn)
                    success = True
                    if self.model_name == 'sugar_ext':
                        self.dic_res['data'][sn_name]['salt2_ext'] = self.source.compute_salt2_params_ext(res_sn['parameters'][2:-2])                 
                except :
                    success = False
                    print ('fit fail for '+sn_name)
                    res_sn, fitted_model_sn = 'fit fail', np.nan
                    self.fit_fail.append(sn_name)
            elif self.model_name == 'salt2' or self.model_name == 'salt2_newerr':
                try :
                    res_sn, fitted_model_sn = self.fit_lc_salt2(self.table_sn)
                    success = True
                except :
                    success = False
                    print ('fit fail for '+sn_name)
                    res_sn, fitted_model_sn = 'fit fail', np.nan
                    self.fit_fail.append(sn_name)
                    
            else :
                raise ValueError('Error model_name have to be salt2 or sugar')

            self.fitted_model.append(fitted_model_sn)
            self.dic_res['data'][sn_name]['data_table'] = np.array(self.table_sn)
            self.dic_res['data'][sn_name]['res'] = res_sn
            self.dic_res['data'][sn_name]['zhl'] = self.zhl
            self.dic_res['data'][sn_name]['mwebv'] = self.mwebv 
            self.dic_res['data'][sn_name]['zerr'] = self.zerr
            self.dic_res['data'][sn_name]['zcmb'] = self.zcmb                
            if self.sample == 'jla':
                self.dic_res['data'][sn_name]['biascor'] = self.biascor
        self.fitting_sample = False
       
    def write_result(self, specific_file=None):
        """
        """
        if specific_file==None:
            if self.dic_res == None:
                raise ValueError('fitted sample needed to write result')
            elif self.sample =='SNf':
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_'+self.strfilters+'_'+self.model_name+'.pkl','wb')
                    pkl.dump(self.dic_res, File)            
            elif self.sample =='csp': 
                if type(self.filter_drop_csp) == str:
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_drop'+self.filter_drop_csp+'_'+self.model_name+'.pkl','wb')
                    pkl.dump(self.dic_res, File)
                elif type(self.filter_drop_csp) == list:
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_drop'+str(self.filter_drop_csp)+'_'+self.model_name+'.pkl','wb')
                    pkl.dump(self.dic_res, File)                
                else: 
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_'+self.model_name+'.pkl','wb')
                    pkl.dump(self.dic_res, File)
            else:
                if self.sub_sample == None:
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_'+self.model_name+'.pkl','wb')
                    pkl.dump(self.dic_res, File)
                elif type(self.sub_sample) == str :
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sub_sample+'_'+self.model_name+'.pkl','wb')
                    pkl.dump(self.dic_res, File)       
                elif type(self.sub_sample) == list :
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'list_'+self.model_name+'.pkl','wb')
                    pkl.dump(self.dic_res, File) 
        else:
            File = open(specific_file, 'wb')
            pkl.dump(self.dic_res, File)
        File.close()


    
