#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:51:58 2019

@author: florian
"""
import sncosmo
import numpy as np
from scipy.interpolate import interp2d
import os 
import abc
from textwrap import dedent
from math import ceil
from iminuit import Minuit
from sncosmo.bandpasses import get_bandpass, Bandpass
from sncosmo.magsystems import get_magsystem
from sncosmo.salt2utils import BicubicInterpolator
from sncosmo.utils import integration_grid
from sncosmo.constants import HC_ERG_AA, MODEL_BANDFLUX_SPACING
from .constant import wl_min_sug
from .load_sugext import register_SUGAR_ext
from .builtins import register_SNf_bands_width, mag_sys_SNF_width,  builtins_jla_bandpasses, mag_sys_jla
from .read_jla import bandpass_interpolators
import pickle as pkl
        
class fit_global_sample_lc(object):
    
    def __init__(self, data, salt2_res_path='../../sugar_analysis_data/extension_training/tlimsug-8salt2-15/resfitlc_snls_salt2.pkl',
                 bands_used=[ 'jla_MEGACAMPSF::g', 'jla_MEGACAMPSF::r','jla_MEGACAMPSF::i', 'jla_MEGACAMPSF::z' ],
                 model_dir = '../../sugar_model/', cov_band=False, radius = None,
                 mod_errfile='../../sugar_analysis_data/err_mod_training/sugar/err_mod_4node.dat', sad_path = '../../',
                 phase_node=[-5., 0., 15., 40.], model='sugar',
                 fit_err_mod=False, fit_A=False, sample='jla',
                 init_A=None, write_chi2=True):
        """
        """
        if model == 'sugar':
            self.model_name = 'sugar'
            self.nb_mod_param = 4
            self.param_mod = ['Xgr', 'q1', 'q3', 'A']
        elif model == 'salt2':
            self.model_name = 'salt2'
            self.nb_mod_param = 3
            self.param_mod = ['x0', 'x1', 'c']
        else:
            raise ValueError('model should be sugar or salt2')
        self.init_A = init_A
        self.write_chi2 = write_chi2
        self.salt2_res_path = salt2_res_path
        self.bands_used = bands_used
        self.phase_node = phase_node
        self.sad_path = sad_path
        register_SUGAR_ext(modeldir=model_dir, mod_errfile=mod_errfile,
                       transmat_path='trans_matrix_init.pkl', version='0.0')
        self.source = sncosmo.get_source('sugar_ext')
        self.fit_err_mod = fit_err_mod
        self.data = data
        self.fit_A = fit_A
        self.cov_band = cov_band
        self.chi_sn = {}
        self.radius = radius
        if self.radius is not None:
            self.radius = pkl.load(open(self.radius, 'rb'))
        if sample == 'jla':
            builtins_jla_bandpasses(sad_path=self.sad_path)
            mag_sys_jla(sad_path=self.sad_path)
            
    def compute_A(self, coeff_A):
        A = np.zeros((5,3))
        for i in range(len(A)):
            A[i,:] = coeff_A[i*3:3*i+3]       
        return A
        
    def err_mod_spline(self, node_errmod):
        self.wave_bin = np.zeros(len(self.bands_used))
        node_array =  np.zeros((len(self.wave_bin), len(self.phase_node)))
        for l in range(len(self.phase_node)):   
            for i in range(len(self.wave_bin)):
                if l ==0:
                    f = sncosmo.get_bandpass(self.bands_used[i])
                    self.wave_bin[i] = float(f.wave_eff)
                    
                node_array[i,l] = node_errmod[l+len(self.phase_node)*i]
        self.spline = interp2d(self.phase_node, self.wave_bin, node_array) 
        
    def intrinsic_dispertion(self, p, wl):
        return self.spline(p, wl)[0]
        
    def chi2(self, param_sug, coeff_A, node_errmod):

        if self.fit_A:
            self.source.update_trans_matrix(self.compute_A(coeff_A))
        dust = sncosmo.CCM89Dust()
        self.model = sncosmo.Model(source=self.source, 
                              effects=[dust], 
                              effect_names=['mw'], 
                              effect_frames=['obs']) 
        if self.fit_err_mod :
            self.err_mod_spline(node_errmod)
        chi = 0
        for j, sn_name in enumerate(self.data.keys()):
            par = param_sug[j*self.nb_mod_param:j*self.nb_mod_param+self.nb_mod_param]
            chi += self.chi_per_sn(sn_name, par)
        return chi

    def get_chi2(self, parameters):
        param_sug = parameters[:len(self.data.keys())*self.nb_mod_param]
        if self.fit_A:
            coeff_A = parameters[len(self.data.keys())*self.nb_mod_param:len(self.data.keys())*self.nb_mod_param+15]        
        else:
            coeff_A = None        
        if self.fit_err_mod:
            node_errmod = parameters[-len(self.phase_node):]
        else:
            node_errmod = None
        chi2 = self.chi2(param_sug, coeff_A, node_errmod)
        if self.write_chi2:
            self.file_chi2.write(str(chi2)+'\n')
        return chi2

    def change_model_params(self, parameters):  
        self.model.set(z=float(self.data[self.sn_name]['res']['parameters'][0]))
        self.model.set(mwebv=self.data[self.sn_name]['mwebv'])
        self.model.set(t0=float(self.data[self.sn_name]['res']['parameters'][1]))

        if self.model_name == 'sugar' :
            self.model.set(Xgr=parameters[0])
            self.model.set(q1=parameters[1])
            self.model.set(q3=parameters[2])
            self.model.set(A=parameters[3])
        elif self.model_name == 'salt2':
            self.model.set(x0=parameters[0])      
            self.model.set(x1=parameters[1])
            self.model.set(c=parameters[2])
            
    def chi_per_sn(self, sn_name, parameters):   
        self.sn_name = sn_name
        self.change_model_params(parameters)
        time_obs = self.data[self.sn_name]['data_table']['time']
        bands = self.data[self.sn_name]['data_table']['band']
        if self.radius is not None:
            rad = self.radius[sn_name]
            for fname in self.bands_used:
                if fname in bands:
                    name = fname[4:]
                    filt = bandpass_interpolators(name, rad)
                    wlen = filt[0]
                    tran = filt[1]
                    fi = sncosmo.Bandpass(wlen, tran, name=fname)
                    sncosmo.registry.register(fi, force=True)
                    bands_ab = {'jla_SDSS::u': ('jla_AB_B12_0',  0.06791),
                            'jla_SDSS::g': ('jla_AB_B12_0', -0.02028),
                            'jla_SDSS::r': ('jla_AB_B12_0', -0.00493),
                            'jla_SDSS::i': ('jla_AB_B12_0', -0.01780),
                            'jla_SDSS::z': ('jla_AB_B12_0', -0.01015),
                            'jla_MEGACAMPSF::u': ('jla_AB_B12_0',  0),
                            'jla_MEGACAMPSF::g': ('jla_AB_B12_0', 0),
                            'jla_MEGACAMPSF::r': ('jla_AB_B12_0', 0),
                            'jla_MEGACAMPSF::i': ('jla_AB_B12_0', 0),
                            'jla_MEGACAMPSF::z': ('jla_AB_B12_0', 0)}
                    sncosmo.registry.register(sncosmo.CompositeMagSystem(bands=bands_ab),'jla_AB_B12', force=True)
        phase = (time_obs- self.data[self.sn_name]['res']['parameters'][1]) / (1 + self.data[self.sn_name]['res']['parameters'][0])  
        data_flux = self.data[self.sn_name]['data_table']['flux']
        data_fluxerr = self.data[self.sn_name]['data_table']['fluxerr']
        zp = self.data[self.sn_name]['data_table']['zp']
        zpsys = self.data[self.sn_name]['data_table']['zpsys']
        
 
        model_flux = self.model.bandflux(bands, time_obs, zp=zp, zpsys=zpsys)
        residuals = data_flux - model_flux
        if self.cov_band:
            cov = self.data[self.sn_name]['data_table']['fluxcov']
            if self.fit_err_mod :           
                weff = np.zeros_like(bands)
                for i, b in enumerate(bands):
                    band = sncosmo.get_bandpass(b) 
                    weff = band.wave_eff/(1 + self.data[self.sn_name]['res']['parameters'][0])
                    cov += np.diag((self.intrinsic_dispertion(phase, weff)*data_flux)**2)
            invcov = np.linalg.pinv(cov)
            chi_sn = np.dot(np.dot(residuals, invcov), residuals)
            self.chi_sn[sn_name] = chi_sn
            return chi_sn   
        else:           
            if self.fit_err_mod :           
                weff = np.zeros_like(bands)
                for i, b in enumerate(bands):
                    band = sncosmo.get_bandpass(b) 
                    weff = band.wave_eff/(1 + self.data[self.sn_name]['res']['parameters'][0])
                w = 1/((data_fluxerr)**2 + (self.intrinsic_dispertion(phase, weff)*model_flux)**2) 
                counter_term = 0 #a remettre 
                log_det_cov = np.sum(np.log(w))
                L = - log_det_cov + chi_sn + counter_term
                self.chi_sn[sn_name] = chi_sn
                return L
            else:
                w = 1/(data_fluxerr)**2
                chi_sn = np.sum(residuals**2*w)
                self.chi_sn[sn_name] = chi_sn
                return chi_sn

    def _init_minimize(self):
        self.param_name = []
        self.init_val = {}
        for sn in self.data.keys():
            for n, par in enumerate(self.param_mod):
                self.param_name.append(sn+par)
                self.init_val[sn+par] = self.data[sn]['res']['parameters'][n+2]      
                self.init_val['error_' + sn+par] = 0.1 * self.data[sn]['res']['parameters'][n+2] 
        if self.fit_A:
            for i in range(15):
                 self.param_name.append('A'+str(i))
                 self.init_val['A'+str(i)] = self.init_A[i]
                 self.init_val['error_' +'A'+str(i)] = 0.1 
        if self.fit_err_mod:
            for j in range(len(self.phase_node)):
                self.param_name.append('node_'+str(j))
                self.init_val['node_'+str(j)] = 0.1
                self.init_val['error_' +'node_'+str(j)] = 0.01 
            
    def minimize(self, file_chi2_path='../../sugar_analysis_data/extension_training/chi2.txt'):
        self._init_minimize()
        if self.write_chi2:
            self.file_chi2 = open(file_chi2_path, 'w')
        self.m = Minuit(self.get_chi2, use_array_call=True, errordef=1.,
                        forced_parameters=self.param_name, **self.init_val)
        self.m.migrad()
        if self.write_chi2:
            self.file_chi2.close()
    
    def write_res(self, 
                  output_path='../../sugar_analysis_data/extension_training/global_fit_SDSS.pkl',
                  chisn_path='../../sugar_analysis_data/extension_training/chipersn.pkl'):
        self.dic_res = {}
        values = self.m.values
        for i, sn in enumerate(self.data.keys()):
            self.dic_res[sn] = {}
            for n, par in enumerate(self.param_mod):
                self.dic_res[sn][par] = values[sn+par]
                self.dic_res[sn]['err_'+par] = self.m.errors[sn+par]
        if self.fit_A:
            coeff_A = []
            for i in range(15):
                coeff_A.append(values['A'+str(i)])
            coeff_A = np.array(coeff_A)
            self.trans_mat_minuit = self.compute_A(coeff_A)        
            self.dic_res['A'] = self.trans_mat_minuit
        if self.fit_err_mod:
            err_mod = []
            for j in range(len(self.phase_node)):
                err_mod.append(values['node_'+str(j)])
            self.dic_res['err_mod'] = np.array(err_mod)
        File = open(output_path, 'wb')
        pkl.dump(self.dic_res, File)
        File.close()
        File_chi_sn = open(chisn_path, 'wb')
        pkl.dump(self.chi_sn, File_chi_sn)
        File_chi_sn.close()
        


        
        
