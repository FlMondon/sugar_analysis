#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 10:33:31 2018

@author: florian
"""

import sncosmo
import builtins_SNF as Build_SNF
import numpy as np
import copy
import constant as cst
import read_SNf_data as RSNf
import cPickle as pkl

SNf_path = '../../sugar_analysis_data/META-CABALLO2.pkl'
output_path = '../../'

class LC_Fitter(object):
    
    def __init__(self, model_name='sugar', sample='SNf', data=None , t0_fix=False, sub_sample=None, modelcov=False):
        
        self.model_name = model_name
        self.t0_fix = t0_fix
        self.dic_res = None
        self.sub_sample = sub_sample
        self.sample = sample
        self.modelcov = modelcov
        Build_SNF.register_SUGAR()
        
        if self.sample=='SNf':
            self.meta = pkl.load(open(SNf_path))
            self.data = {}
            if sub_sample==None:
                for sn_name in self.meta.keys():
                    if self.meta[sn_name]['idr.subset'] == 'training' or self.meta[sn_name]['idr.subset'] == 'validation':
                        self.data[sn_name] = self.meta[sn_name]
            else :
                for sn_name in self.sub_sample:
                    self.data[sn_name] = self.meta[sn_name]
                    
            self.filters = ['BSNf','VSNf','RSNf']
            self.errorscale = True
            self.width= 10
            Build_SNF.register_SNf_bands_width(width=self.width)
            Build_SNF.mag_sys_SNF_width(width=self.width)
            
        elif sample=='jla':
            raise ValueError('jla not implemented yet')
        else:
            Build_SNF.register_SNf_bands_width(width=self.width)
            Build_SNF.mag_sys_SNF_width(width=self.width)
            print 'Warning : Mag sys and band used in your data have to be register in sncosmo our in builtins_SNF'
            try:
                self.data = pkl.load(open(data))
                if type(self.data) is not dict:
                    raise ValueError('data type have to be a dictionary')
            except:
                raise ValueError('data type have to be a dictionary')
    
    def sn_data(self, sn_name):
       
        if self.sample == 'SNf':
            self.table_sn, self.zhl, self.zerr, self.zcmb, self.mwebv, self.daymax = RSNf.read_meta_SNF(self.data, sn_name,
                                               filters=self.filters, 
                                               model=self.model_name, 
                                               errorscale=self.errorscale,
                                               width=self.width)
        elif self.sample == 'jla':
            raise ValueError('jla not implemented yet')
            
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
            
    
    def fit_lc_sugar(self, data, zhl=None,  mwebv=None):
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
            
        self.model.set(mwebv=self.mwebv)
        self.model.set(z=self.zhl)
        
        print 'initialisation'            

        self.model.set(q1=0.0)
        self.model.set(q2=0.0)
        self.model.set(q3=0.0)     
        self.model.set(Xgr=1.0e-15)
        
        if self.t0_fix:
            self.model.set(t0=self.daymax)
            res, fitted_model = sncosmo.fit_lc(data, 
                                               self.model, 
                                               ['A','Xgr'], 
                                               modelcov=False)
            res, fitted_model = sncosmo.fit_lc(data, 
                                               self.model, 
                                               ['q1', 
                                                'q2', 
                                                'q3', 
                                                'A', 
                                                'Xgr'], 
                                               modelcov  = self.modelcov)
            print res.chisq
        else:
            res, fitted_model = sncosmo.fit_lc(data, 
                                               self.model, 
                                               ['t0','A','Xgr'], 
                                               modelcov=False)
    
            print res.parameters[1] 
            print 'first iteration'
            chi2 = res.chisq
            chi2p = chi2*2
            m=0
            while chi2 < chi2p and m < 10:
    
                print m
                if m > 0:
                    resp = res
                    fitted_modelp = fitted_model
                m += 1
                t_peak = fitted_model.parameters[1]
                t1 = t_peak + cst.t_min_sug*(1 + self.model.get('z'))
                t2 = t_peak + cst.t_max_sug*(1 + self.model.get('z'))
                            
                A=[]
                data_new = copy.deepcopy(data)
                for i in range(len(data_new)):                    
                    if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                        A.append(i)
                A=np.array(A)
                for i in range(len(A)):
                    data_new.remove_row(A[i])
                    A-=1   
    
    
                print "WARNING: Sugar don't have model covariance for the moment"
                
                self.model.set(q1=res.parameters[3])
                self.model.set(q2=res.parameters[4])
                self.model.set(q3=res.parameters[5])
                self.model.set(A=res.parameters[6])
                self.model.set(Xgr=res.parameters[2])
                self.model.set(t0=res.parameters[1])                        
                res, fitted_model = sncosmo.fit_lc(data_new, 
                                                   self.model, 
                                                   ['t0', 
                                                    'q1', 
                                                    'q2', 
                                                    'q3', 
                                                    'A', 
                                                    'Xgr'], 
                                                   modelcov  = self.modelcov,
                                                   bounds = {'t0' : [res.parameters[1]-30,res.parameters[1]+30]})
                
                chi2p = chi2
                chi2 = res.chisq
                print chi2p, chi2
            #final results
            res = resp
            fitted_model = fitted_modelp
        return res, fitted_model    

    def fit_lc_salt2(self, data, zhl=None,  mwebv=None ):
        """
        """
        if not self.fitting_sample :
            self.source = sncosmo.get_source('salt2')
            dust = sncosmo.CCM89Dust()
            self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs']) 
            self.mwebv = mwebv
            self.zhl = zhl
            
        self.model.set(mwebv=self.mwebv)
        self.model.set(z=self.zhl)
        print 'initialisation'            
        res, fitted_model = sncosmo.fit_lc(data, 
                                           self.model, 
                                           ['t0','c','x0'], 
                                           modelcov=False)

        print res.parameters[1] 
        print 'first iteration'
        chi2 = res.chisq
        print chi2
        chi2p = chi2*2
        m=0
    
        while chi2 < chi2p and m < 10:

            print m
            if m > 0:
                resp = res
                fitted_modelp = fitted_model
            m += 1
            t_peak = fitted_model.parameters[1]
            #print t_peak,fitted_model.parameters[4]

            t1 = t_peak + cst.t_min_sug*(1 + self.model.get('z'))
            t2 = t_peak + cst.t_max_sug*(1 + self.model.get('z'))
                        
            A=[]
            data_new = copy.deepcopy(data)
            for i in range(len(data_new)):                    
                if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                    A.append(i)
            A=np.array(A)
            for i in range(len(A)):
                data_new.remove_row(A[i])
                A-=1   


            print "WARNING: Sugar don't have model covariance for the moment"
            modelcov = False
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
                                               modelcov  = modelcov,
                                               bounds = {'t0' : [res.parameters[1]-30,res.parameters[1]+30]})
            
            chi2p = chi2
            chi2 = res.chisq
            print chi2p, chi2
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
        self.fit_fail = []
        self.fitted_model = []
        self.fitting_sample = True
        self.source = sncosmo.get_source(self.model_name)
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs'])
        for sn_name in self.data.keys():
            self.sn_data(sn_name)
            print sn_name
            if self.model_name == 'sugar':
                try:
                    res_sn, fitted_model_sn = self.fit_lc_sugar(self.table_sn)
                except :
                    print 'fit fail for '+sn_name
                    res_sn, fitted_model_sn = 'fit fail', np.nan
                    self.fit_fail.append(sn_name)
            elif self.model_name == 'salt2':
                try :
                    res_sn, fitted_model_sn = self.fit_lc_salt2(self.table_sn)
                except :
                    print 'fit fail for '+sn_name
                    res_sn, fitted_model_sn = 'fit fail', np.nan
                    self.fit_fail.append(sn_name)
            else :
                raise ValueError('Error model_name have to be salt2 or sugar')
            self.fitted_model.append(fitted_model_sn)
            self.dic_res['data'][sn_name] = {}
            self.dic_res['data'][sn_name]['res'] = res_sn
            self.dic_res['data'][sn_name]['zhl'] = self.zhl
            self.dic_res['data'][sn_name]['zerr'] = self.zerr
            self.dic_res['data'][sn_name]['zcmb'] = self.zcmb
            self.dic_res['data'][sn_name]['mwebv'] = self.mwebv 
        self.fitting_sample = False
       
    def write_result(self):
        """
        """
        if self.dic_res == None:
            raise ValueError('fitted sample needed to write result')
        else:
            File = open(output_path+'sugar_analysis_data/resfitlc_'+self.sample+'_'+self.model_name+'.pkl','w')
            
            pkl.dump(self.dic_res, File)
        

    
    
