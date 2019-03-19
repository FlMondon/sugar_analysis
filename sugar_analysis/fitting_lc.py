#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 10:33:31 2018

@author: florian
"""

import sncosmo
import copy
import numpy as np
import os
try:
   import cPickle as pkl
except ModuleNotFoundError:
   import pickle as pkl
   
from .builtins import register_SNf_bands_width, mag_sys_SNF_width,  builtins_jla_bandpasses, mag_sys_jla
from .load_sugar import register_SUGAR
from .read_jla import read_lc_jla
from .cosmo_tools import distance_modulus_th
from .constant import t_min_sug, t_max_sug, t_min_salt2, t_max_salt2
from .data_table import build_data
from .data_table import read_csp


   
class LC_Fitter(object):
    
    def __init__(self, model_name='sugar', sample='SNf', data=None, 
                 t0_fix=False, sub_sample=None, modelcov=False,
                 filters=['new_fU_10','fB_10','fV_10','fR_10','new_fI_10'], 
                 filter_drop_csp = None, sad_path = '../../',
                 modeldir='../../sugar_model/',
                 width=10, param_sug_fix =[]):
        
        
        self.fitting_sample =  False
        self.model_name = model_name
        self.t0_fix = t0_fix
        self.sad_path = sad_path
        self.filter_drop_csp = filter_drop_csp
        self.dic_res = None
        self.sub_sample = sub_sample
        self.sample = sample
        self.modeldir = modeldir
        self.modelcov = modelcov
        self.filters = filters
        self.strfilters = str(len(filters))+filters[0]+filters[len(filters)-1]
        register_SUGAR(modeldir=self.modeldir)
        self.param_sug = ['q1', 'q2', 'q3']
        for psf in param_sug_fix:
            self.param_sug.remove(psf)

                
        
        if self.sample=='SNf':
            self.bd = build_data(sad_path=self.sad_path, modeldir=self.modeldir)
            self.meta = pkl.load(open(self.sad_path+'sugar_analysis_data/META-CABALLO2.pkl'))
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
            register_SNf_bands_width(width=self.width)
            mag_sys_SNF_width(width=self.width)
            
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
            builtins_jla_bandpasses()
            mag_sys_jla()
            self.rcsp = read_csp()
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
            self.source = sncosmo.get_source('sugar')
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

        print ('initialisation')            

        self.model.set(q1=0.0)
        self.model.set(q2=0.0)
        self.model.set(q3=0.0)
        self.model.set(A=0.1)
        self.model.set(Xgr=Xgr_init)
        if self.t0_fix:
            self.model.set(t0=self.daymax)
            res, fitted_model = sncosmo.fit_lc(data, 
                                               self.model, 
                                               self.param_sug+['A', 'Xgr'], 
                                               modelcov  = self.modelcov)
            print (res.chisq)
        else:
            res, fitted_model = sncosmo.fit_lc(data, 
                                               self.model, 
                                               ['t0','A','Xgr'], 
                                               modelcov=False)
    
            chi2 = res.chisq
            chi2p = chi2*2
            m=0
            recovery = True
            t0_init = res.parameters[1] 
            Xgr_init = res.parameters[2]
            self.model.set(A=res.parameters[6])
            self.model.set(Xgr=res.parameters[2])
            self.model.set(t0=res.parameters[1]) 
            res, fitted_model = sncosmo.fit_lc(data, 
                                               self.model, 
                                               self.param_sug, 
                                               modelcov=False)
            
            print ('first iteration')
            while chi2 < chi2p and m < 10:
    
                print (m)
                if m > 0:
                    resp = res
                    fitted_modelp = fitted_model
                t_peak = fitted_model.parameters[1]
                t1 = t_peak + t_min_sug*(1 + self.model.get('z'))
                t2 = t_peak + t_max_sug*(1 + self.model.get('z'))
                            
                A=[]
                data_new = copy.deepcopy(data)
                for i in range(len(data_new)):                    
                    if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                        A.append(i)
                A=np.array(A)
                for i in range(len(A)):
                    data_new.remove_row(A[i])
                    A-=1   
    
    
                Warning("Sugar don't have model covariance for the moment")
                
                self.model.set(q1=res.parameters[3])
                self.model.set(q2=res.parameters[4])
                self.model.set(q3=res.parameters[5])
                self.model.set(A=res.parameters[6])
                self.model.set(Xgr=res.parameters[2])
                self.model.set(t0=res.parameters[1])                        
                res, fitted_model = sncosmo.fit_lc(data_new, 
                                                   self.model, 
                                                   ['t0']+self.param_sug+['A', 'Xgr'], 
                                                   modelcov  = self.modelcov,
                                                   bounds = {'t0' : [res.parameters[1]-30,res.parameters[1]+30]})
            
                chi2p = chi2
                chi2 = res.chisq
                print (chi2p, chi2)
                if m == 0 and chi2 > chi2p :
                    if recovery:
                        print ('First iteration fail try to recovery :')
                        chi2 = self.recovery_sugar(data, t0_init, Xgr_init, res, chi2p)
                        recovery = False 
                    else:
                        raise ValueError('Error initial value better than first iteration')
                else:
                    m += 1
            
            #final results
            res = resp
            fitted_model = fitted_modelp
        return res, fitted_model 
    
    def recovery_sugar(self, data, t0_init, Xgr_init, res, chi2p):
        n=0
        t0_iter=t0_init-5
        while res.chisq > chi2p and n < 10:
            print (n)
            n+=1
    
            self.model.set(Xgr=Xgr_init)
            self.model.set(q1=0.0)
            self.model.set(q2=0.0)
            self.model.set(q3=0.0)
            self.model.set(A=0.2)
            self.model.set(t0=t0_iter)   
            try:
                res, fitted_model = sncosmo.fit_lc(data, 
                                           self.model, 
                                           ['t0']+self.param_sug+['A', 'Xgr'], 
                                           modelcov  = self.modelcov,
                                           bounds = {'t0' : [res.parameters[1]-5,res.parameters[1]+5]}) 
            except:
                'Fit fail in recovery iteration %f! Next iteration'%n
            t0_iter += 1
            print (res.chisq)
        if n == 10:
            raise ValueError('recovery fail')
        else:
            print ('successful recovery retry to fit')
        
        return res.chisq
            
    def fit_lc_salt2(self, data, zhl=None,  mwebv=None ):
        """
        """
        if not self.fitting_sample :
            self.source = sncosmo.get_source('salt2', version='2.4')
            dust = sncosmo.CCM89Dust()
            self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs'])         
        self.model.set(mwebv=self.mwebv)
        self.model.set(z=self.zhl)
        print ('initialisation')            
        res, fitted_model = sncosmo.fit_lc(data, 
                                           self.model, 
                                           ['t0','c','x0'], 
                                           modelcov=False)

        print (res.parameters[1])
        print ('first iteration')
        chi2 = res.chisq
        print (chi2)
        chi2p = chi2*2
        m=0
    
    
        while chi2 < chi2p and m < 20:

            print (m)
            if m > 0:
                resp = res
                fitted_modelp = fitted_model
            m += 1
            t_peak = fitted_model.parameters[1]
            #print t_peak,fitted_model.parameters[4]

            t1 = t_peak + t_min_salt2*(1 + self.model.get('z'))
            t2 = t_peak + t_max_salt2*(1 + self.model.get('z'))
                        
            A=[]
            data_new = copy.deepcopy(data)
            for i in range(len(data_new)):                    
                if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                    A.append(i)
            A=np.array(A)
            for i in range(len(A)):
                data_new.remove_row(A[i])
                A-=1   

            self.model.set(x1=res.parameters[3])
            self.model.set(c=res.parameters[4])
            self.model.set(x0=res.parameters[2])
            self.model.set(t0=res.parameters[1])                        
            res, fitted_model = sncosmo.fit_lc(data_new, 
                                               self.model, 
                                               ['t0', 
                                                'x1',
                                                'c', 
                                                'x0'], 
                                               modelcov  = True,
                                               bounds = {'t0' : [res.parameters[1]-30,res.parameters[1]+30]})
            
            chi2p = chi2
            chi2 = res.chisq
            print (chi2p, chi2)
        #final results
        res = resp
        fitted_model = fitted_modelp
        return res, fitted_model
    
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
        if self.model_name== 'salt2':
            self.source = sncosmo.get_source(self.model_name, version='2.4')
        else:
            self.source = sncosmo.get_source(self.model_name)
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs'])
        for sn_name in self.data:
            self.sn_data(sn_name)
            print (sn_name)
            if self.model_name == 'sugar':
                try:
                    res_sn, fitted_model_sn = self.fit_lc_sugar(self.table_sn)
                    success = True
                except :
                    success = False
                    print ('fit fail for '+sn_name)
                    res_sn, fitted_model_sn = 'fit fail', np.nan
                    self.fit_fail.append(sn_name)
            elif self.model_name == 'salt2':
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
            self.dic_res['data'][sn_name] = {}
#            self.dic_res['data'][sn_name]['fitted_model'] = fitted_model_sn
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
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_'+self.strfilters+'_'+self.model_name+'.pkl','w')
                    pkl.dump(self.dic_res, File)            
            elif self.sample =='csp': 
                if type(self.filter_drop_csp) == str:
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_drop'+self.filter_drop_csp+'_'+self.model_name+'.pkl','w')
                    pkl.dump(self.dic_res, File)
                elif type(self.filter_drop_csp) == list:
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_drop'+str(self.filter_drop_csp)+'_'+self.model_name+'.pkl','w')
                    pkl.dump(self.dic_res, File)                
                else: 
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_'+self.model_name+'.pkl','w')
                    pkl.dump(self.dic_res, File)
            else:
                if self.sub_sample == None:
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'_'+self.model_name+'.pkl','w')
                    pkl.dump(self.dic_res, File)
                elif type(self.sub_sample) == str :
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sub_sample+'_'+self.model_name+'.pkl','w')
                    pkl.dump(self.dic_res, File)       
                elif type(self.sub_sample) == list :
                    File = open(self.sad_path+'sugar_analysis_data/resfitlc_'+self.sample+'list_'+self.model_name+'.pkl','w')
                    pkl.dump(self.dic_res, File) 
        else:
            pkl.dump(self.dic_res, open(specific_file, 'w'))


    
