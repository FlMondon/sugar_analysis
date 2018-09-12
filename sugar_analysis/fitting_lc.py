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
output_path = '../../sugar_analysis_data/'

class LC_Fitter(object):
    
    def __init__(self, model_name='sugar', sample='SNf', data=None , t0_fix=False, sub_sample=None):
        
        self.model_name = model_name
        self.t0_fix = t0_fix
        self.dic_res = None
        self.sub_sample = sub_sample
        self.sample = sample
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
            self.table_sn, self.zhl, self.zcmb, self.mwebv, self.daymax = RSNf.read_meta_SNF(self.data, sn_name,
                                               filters=self.filters, 
                                               model=self.model_name, 
                                               errorscale=self.errorscale,
                                               width=self.width)
        elif self.sample == 'jla':
            raise ValueError('jla not implemented yet')
            
        else:
            self.table_sn = self.data[sn_name]['table']
            self.zhl = self.data[sn_name]['zhl']
            self.zcmb = self.data[sn_name]['zcmb']
            self.mwebv = self.data[sn_name]['mwebv']
            if self.t0_fix== True:
                try : 
                    self.daymax = self.data[sn_name]['daymax']
                except:
                    raise ValueError('t0_fix is True: need daymax in the dictionary data')
            
    
    def fit_lc_sugar(self, data ):
        """
        """
        if not self.fitting_sample :
            self.source = sncosmo.get_source('sugar')
            dust = sncosmo.CCM89Dust()
            self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs']) 

        self.model.set(mwebv=self.mwebv)
        self.model.set(z=self.zhl)
        
        print 'initialisation'            

        self.model.set(q1=0.0)
        self.model.set(q2=0.0)
        self.model.set(q3=0.0)     
        self.model.set(Xgr=1.0e-15)
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
            modelcov = False
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
                                               modelcov  = modelcov,
                                               bounds = {'t0' : [res.parameters[1]-30,res.parameters[1]+30]})
            
            chi2p = chi2
            chi2 = res.chisq
            print chi2p, chi2
        #final results
        res = resp
        fitted_model = fitted_modelp
        return res, fitted_model    

    def fit_lc_salt2(self, data ):
        """
        """
        if not self.fitting_sample :
            self.source = sncosmo.get_source('salt2')
            dust = sncosmo.CCM89Dust()
            self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs']) 
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
                                               modelcov  = modelcov)
            
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
        self.fit_fail = []
        self.fitting_sample = True
        self.source = sncosmo.get_source(self.model_name)
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs'])
        for sn_name in self.data.keys():
            self.sn_data(sn_name)
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
            self.dic_res[sn_name] = {}
            self.dic_res[sn_name]['res'] = res_sn
            self.dic_res[sn_name]['fitted_model'] = fitted_model_sn
        self.fitting_sample = False
       
    def write_result(self):
        """
        """
        if self.dic_res == None:
            raise ValueError('fitted sample needed to write result')
        else:
            pkl.dump(self.dic_res, output_path+'resfitlc_'+self.sample+'_'+self.model_name)
        
        
        
        
        
    def mB_determination(res):
        """
        determine mb restframe 
        """
        scale_factor = 10**-12
        
        #interpolation of TB and Trest
        filt2 = np.genfromtxt(jla_path+'Instruments/SNLS3-Landolt-model/sb-shifted.dat')
        wlen = filt2[:,0]
        tran = filt2[:,1]
        splB = Spline1d(wlen, tran, k=1,ext = 1)
    
    
    
        #interpolation of ref spectrum
        data = np.genfromtxt(jla_path+'MagSys/bd_17d4708_stisnic_002.ascii')
        dispersion = data[:,0]
        flux_density = data[:,1]
    #    fits = pyfits.open(sugar_analysis_data+'Vega.fits')
    #    fit = fits[1]
    #    wl = np.zeros(len(fit.data))
    #    flux = np.zeros(len(fit.data))
    #    for i in range(len(fit.data)):
    #         wl[i] = fit.data[i][0]
    #         flux[i] = fit.data[i][1]
    
        splref = Spline1d(dispersion, flux_density, k=1,ext = 1)
    
      
        #interpolation of the spectrum flux
        template_0 = np.genfromtxt(jla_path+'salt2-4/salt2_template_0.dat')    
        template_1 = np.genfromtxt(jla_path+'salt2-4/salt2_template_1.dat')
    #    salt2source=sncosmo.SALT2Source('/users/divers/lsst/mondon/hubblefit/sncosmo_jla/salt2-4')
        
        wlM0 = []
        M0 = []
        for i in range(len(template_0[:,0])):
            if template_0[:,0][i] == 0.0:
                 wlM0.append(template_0[:,1][i]) 
                 M0.append(template_0[:,2][i])
        splM0 = Spline1d(wlM0, M0, k=1,ext = 1)
    
        wlM1 = []
        M1 = []
        for i in range(len(template_1[:,0])):
            if template_1[:,0][i] == 0.0:
                wlM1.append(template_1[:,1][i]) 
                M1.append(template_1[:,2][i])
        splM1 = Spline1d(wlM1, M1, k=1,ext = 1)
    
        #computation of the integral
        dt = 100000
        xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
        dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))
    
    #    I1=np.sum((splM0(xs)*10**-12+res.parameters[3]*splM1(xs)*10**-12)*(10**(-0.4*salt2source.colorlaw(xs)*res.parameters[4]))*xs*splB(xs)*dxs)
        I1 = np.sum((splM0(xs)*scale_factor + res.parameters[3]*splM1(xs)*scale_factor)*(10**(-0.4*color_law_salt2(xs)*res.parameters[4]))*xs*splB(xs)*dxs)    
        I2 = np.sum(splref(xs)*xs*splB(xs)*dxs)
    #    print I1, I2   
        
        
        
        #computation of mb
        mref = 9.907
        mb = -2.5*np.log10(res.parameters[2]*(I1/I2))+mref
    
    
    def color_law_salt2(self, wl):
        B_wl = 4302.57
        V_wl = 5428.55
        l = (wl-B_wl)/(V_wl-B_wl)
        l_lo = (2800.-B_wl)/(V_wl-B_wl)
        l_hi = (7000.-B_wl)/(V_wl-B_wl)
        a = -0.504294
        b = 0.787691
        c = -0.461715
        d = 0.0815619
        cst = 1-(a+b+c+d)
        cl = []
        for i in range (len(l)):
            if l[i] > l_hi:
                cl.append(-(cst*l_hi+l_hi**2*a+l_hi**3*b+l_hi**4*c+l_hi**5*d+(cst+2*l_hi*a+3*l_hi**2*b+4*l_hi**3*c+5*l_hi**4*d)*(l[i]-l_hi)))
            if l[i] < l_lo:
                cl.append(-(cst*l_lo+l_lo**2*a+l_lo**3*b+l_lo**4*c+l_lo**5*d+(cst+2*l_lo*a+3*l_lo**2*b+4*l_lo**3*c+5*l_lo**4*d)*(l[i]-l_lo))) 
            if l[i]>= l_lo and l[i]<= l_hi:
                cl.append(-(cst*l[i]+l[i]**2*a+l[i]**3*b+l[i]**4*c+l[i]**5*d)) 
        return np.array(cl)