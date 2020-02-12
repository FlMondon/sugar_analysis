#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 17:24:17 2019

@author: florian
"""

import sncosmo
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib import pyplot as plt
try:
   import cPickle as pkl
except ModuleNotFoundError:
   import pickle as pkl
from .builtins import register_SNf_bands_width, mag_sys_SNF_width,  builtins_jla_bandpasses, mag_sys_jla
from .load_sugar import register_SUGAR
from .read_jla import read_lc_jla
from .load_sugext import register_SUGAR_ext

import numpy as np
import copy
from astropy.stats import median_absolute_deviation as mad_astropy

def biweight_M(sample,CSTD=6.):
    """
    average using biweight (beers 90)
    """
    M = np.median(sample)
    iterate = [copy.deepcopy(M)]
    mu = (sample-M)/(CSTD*mad_astropy(sample))
    Filtre = (abs(mu)<1)
    up = (sample-M)*((1.-mu**2)**2)
    down = (1.-mu**2)**2
    M = M + np.sum(up[Filtre])/np.sum(down[Filtre])
    iterate.append(copy.deepcopy(M))
    i=1
    while abs((iterate[i-1]-iterate[i])/iterate[i])<0.001:
        mu = (sample-M)/(CSTD*mad_astropy(sample))
        Filtre = (abs(mu)<1)
        up=(sample-M)*((1.-mu**2)**2)
        down = (1.-mu**2)**2
        M = M + np.sum(up[Filtre])/np.sum(down[Filtre])
        iterate.append(copy.deepcopy(M))
        i+=1
        if i == 100 :
            print('voila voila ')
            break
    return M

def biweight_S(sample,CSTD=9.):
    """
    std from biweight
    """
    M = biweight_M(sample)
    mu = (sample-M)/(CSTD*mad_astropy(sample))
    Filtre = (abs(mu)<1)
    up = ((sample-M)**2)*((1.-mu**2)**4)
    down = (1.-mu**2)*(1.-5.*mu**2)
    std = np.sqrt(len(sample))*(np.sqrt(np.sum(up[Filtre]))/abs(np.sum(down[Filtre])))
    return std

class sugar_analysis_plot_lc(object):

    def __init__(self,
                 path_list = ['../../sugar_analysis_data/resfitlc_SNf_5new_fU_10new_fI_10_sugar.pkl',
                                    '../../sugar_analysis_data/resfitlc_SNf_5new_fU_10new_fI_10_salt2.pkl'],
                 sad_path='../../',
                 modeldir='../../sugar_model/', 
                 mod_errfile='../../sugar_analysis_data/err_mod_training/sugar/err_mod_4node.dat',
                 label = ['sugar', 'salt2'],
                 sample='jla', transmat_path='trans_matrix_init.pkl',
                 sys_used = 'jla_AB_B12', version='0.0'): 
        """
        doc a faire 
        """
        if len(path_list) != len(label):
            raise ValueError('path_list and label should have the same dimension')
        
        self.dic_res = {}
        for j, path in enumerate(path_list) :
            print(path)
            dic =  pkl.load(open(path, 'rb'))
            dic_res  = dic['data']
            dic_del = {}
            fit_exclude = []
            for sn_name in dic_res.keys():
                if dic_res[sn_name]['res'] != 'fit fail':
                    dic_del[sn_name] = dic_res[sn_name]
                else:
                    fit_exclude.append(sn_name)
            self.dic_res[label[j]] = dic_del 
        register_SUGAR_ext(modeldir=modeldir, mod_errfile=mod_errfile,
                           transmat_path=transmat_path, version=version)
        if sample == 'jla':
            builtins_jla_bandpasses()
            mag_sys_jla()
        elif sample == 'SNf':
            register_SNf_bands_width(width=10)
            mag_sys_SNF_width(width=10)
        self.sys_used = sys_used
        self.sample = sample
        self.sys = sncosmo.get_magsystem(self.sys_used)
        self.key_label = label
        self.sad_path = sad_path
    
    def build_model(self, params, model_used):
            source = sncosmo.get_source(model_used)
            dust = sncosmo.CCM89Dust()
            model = sncosmo.Model(source=source, 
                      effects=[dust], 
                      effect_names=['mw'], 
                      effect_frames=['obs']) 
            
            model.set(z=float(params[0]))
            model.set(t0 = params[1])
            if model_used == 'sugar':
                model.set(Xgr = float(params[2]))
                model.set(q1 = float(params[3]))
                model.set(q2 = float(params[4]))
                model.set(q3 = float(params[5]))
                model.set(A = float(params[6]))
                model.set(mwebv = float(params[7]))
            elif model_used == 'salt2':
                model.set(x0 = float(params[2]))
                model.set(x1 = float(params[3]))
                model.set(c = float(params[4]))
                model.set(mwebv = float(params[5]))
            elif model_used == 'sugar_ext':
                model.set(Xgr = float(params[2]))
                model.set(q1 = float(params[3]))
                model.set(q3 = float(params[4]))
                model.set(A = float(params[5]))
                model.set(mwebv = float(params[6]))
            return model 
        
    def compute_bin_residuals(self, filters, phase,
                    key_label='sugar', model_used='sugar',
                    total_sn = False, res_path=None):
        """
        compute light curves residuals by filters and by phase
        """
        
        if type(filters) is not list :
            raise ValueError('filters should be a list')
        results = {}
        if total_sn:
            sn = 'total'
            results[sn] = {}
            for f in filters:
                results[sn][f] =  {'residue':[], 'data_err':[],'mag':[], 
                       'time':[], 'zp':[], 'zpsys':[], 'res_bin_phase':{},
                           'var_bin_phase':{}}
                for i in range(len(phase)-1):
                    results[sn][f]['res_bin_phase'][str(phase[i])+'_'+str(phase[i+1])] = []
                    results[sn][f]['var_bin_phase'][str(phase[i])+'_'+str(phase[i+1])] = []

                    
        for sn_name in self.dic_res[key_label].keys():
            if not total_sn:
                sn = sn_name
                results[sn] = {}
                
                for f in filters:
                    results[sn][f] =  {'residue':[], 'data_err':[], 'mag':[],
                           'time':[],'zp':[], 'zpsys':[],'res_bin_phase':{},
                           'var_bin_phase':{}}
                    for i in range(len(phase)-1):
                        results[sn][f]['res_bin_phase'][str(phase[i])+'_'+str(phase[i+1])] = []
                        results[sn][f]['var_bin_phase'][str(phase[i])+'_'+str(phase[i+1])] = []

            params = self.dic_res[key_label][sn_name]['res']['parameters']
            
            if self.sample == 'jla' :
                head, table = read_lc_jla(sn_name, sad_path=self.sad_path, model=model_used)
                self.sys = sncosmo.get_magsystem(self.sys_used)

            else :
                table = self.dic_res[key_label][sn_name]['data_table']

            tmax = float(params[1])      
            model = self.build_model(params, model_used)
            for i, band in enumerate(table['band']):
                band = str(band)
                if band in filters or 'cspv' in band:
                    time = table['time'][i]
                    if time >= tmax + phase[0] and time <= tmax + phase[-1] :
                        try:
                            if self.sample != 'SNf':
                                model_flux = model.bandflux(band, time, zp=table['zp'][i], zpsys=table['zpsys'][i] )
                                model_mag = self.sys.band_flux_to_mag(model_flux, band)
                            else:
                                model_mag = model.bandmag(band, self.sys_used, time)
                            mag = self.sys.band_flux_to_mag(table['flux'][i], band)
                        except:
                            continue
                        if 'cspv' in band:
                            band ='cspv'
                        residuals = mag - model_mag
                        err_mag = table['fluxerr'][i] /table['flux'][i] *1.0857362047581294
                        for i in range(len(phase)-1):
                            if (time-tmax) >= phase[i] and (time-tmax)<=phase[i+1]:
                                results[sn][band]['res_bin_phase'][str(phase[i])+'_'+str(phase[i+1])].append(residuals)
                                results[sn][band]['var_bin_phase'][str(phase[i])+'_'+str(phase[i+1])].append(np.float(err_mag)**2)
        
        if res_path != None:
            File = open(res_path, 'wb')
            pkl.dump(results, File)
        return results
    
    def compute_residuals(self, tlim= [-8, 48], key_label='sugar', 
                          model_used='sugar', res_path=None):
        
        results = {}   
        for sn_name in self.dic_res[key_label].keys():
            results[sn_name] = {}
            params = self.dic_res[key_label][sn_name]['res']['parameters']
            if self.sample == 'jla' :
                head, table = read_lc_jla(sn_name, sad_path=self.sad_path, model=model_used)
                self.sys = sncosmo.get_magsystem(self.sys_used)

            else :
                table = self.dic_res[key_label][sn_name]['data_table']
        
            tmax = float(params[1])      
            model = self.build_model(params, model_used)
            for i, band in enumerate(table['band']):
                band = str(band)
                sncosmo_band = sncosmo.get_bandpass(band)
                if (sncosmo_band.minwave() < model.minwave() or sncosmo_band.maxwave() > model.maxwave()):
                    continue
                if band not in results[sn_name].keys():
                    results[sn_name][band] =  {'residue_mag':[],
                           'residue_flux':[], 'data_err_mag':[], 'mag':[],
                           'data_err_flux':[],  'flux':[],
                           'time':[], 'phase':[], 'zp':[], 'zpsys':[]}                    
                time = table['time'][i]
                if time >= tmax + tlim[0] and time <= tmax + tlim[1] :
                    model_flux = model.bandflux(band, time, zp=table['zp'][i], zpsys=table['zpsys'][i] )
                    residuals_flux = table['flux'][i] - model_flux
                    if model_flux > 0 :
                        model_mag = self.sys.band_flux_to_mag(model_flux, band)
                    else : 
                        model_mag = np.nan
                    if table['flux'][i] > 0 :    
                        mag = self.sys.band_flux_to_mag(table['flux'][i], band)
                    else : 
                        mag = np.nan
                    residuals_mag = mag - model_mag
                else :
                    residuals_mag = np.nan
                    residuals_flux = np.nan
                    if table['flux'][i] > 0.:
                        mag = self.sys.band_flux_to_mag(table['flux'][i], band)
                    else:
                        mag = np.nan
                results[sn_name][band]['phase'].append(time-tmax)
                results[sn_name][band]['time'].append(time)
                residuals_flux = residuals_flux
                results[sn_name][band]['mag'].append(mag)
                results[sn_name][band]['flux'].append(table['flux'][i])
                results[sn_name][band]['residue_mag'].append(residuals_mag)
                results[sn_name][band]['residue_flux'].append(residuals_flux)
                err_mag = table['fluxerr'][i] /table['flux'][i] *1.0857362047581294
                results[sn_name][band]['data_err_mag'].append(err_mag)
                results[sn_name][band]['data_err_flux'].append(table['fluxerr'][i])
                results[sn_name][band]['zp'].append(table['zp'][i])
                results[sn_name][band]['zpsys'].append(table['zpsys'][i])        
        if res_path != None:
                File = open(res_path, 'wb')
                pkl.dump(results, File)
        return results
    
    def prepare_plot(self, font_size=5., 
                     font_family= 'normal',
                     linewidth_axes=1,
                     linewidth_legend=10,
                     xtick_major_size=1,
                     xtick_major_width=1,
                     ytick_major_size=1,
                     ytick_major_width=1,
                     xtick_minor_size=1,
                     xtick_minor_width=1,
                     ytick_minor_size=1,
                     ytick_minor_width=1):
        from matplotlib import rc, rcParams
        rcParams['font.size'] = font_size
        rc('axes', linewidth=linewidth_axes)
        rc('legend', fontsize=linewidth_legend)
        rc('xtick.major', size=xtick_major_size, width=xtick_major_width)
        rc('ytick.major', size=ytick_major_size, width=ytick_major_width)
        rc('xtick.minor', size=xtick_minor_size, width=xtick_minor_width)
        rc('ytick.minor', size=ytick_minor_size, width=ytick_minor_width)

##################casser#####################
#    def plot_distrib(self, filters, phase,  figsize=(15, 25),
#                    BINNING = [np.linspace(-0.5,0.5, 50), 
#                               np.linspace(-0.5,0.5, 50),
#                               np.linspace(-0.5,0.5, 50),
#                               np.linspace(-0.5,0.5, 50),
#                               np.linspace(-0.5,0.5, 50),
#                               np.linspace(-0.5,0.5, 50),
#                               np.linspace(-0.5,0.5, 50)],
#                    NORM=True,
#                    color=['b', 'r', 'black', 'g'],
#                    txt_pos_x=[-0.45, 0.2],
#                    txt_pos_y=[2.5, 1.2],
#                    model_used=None,
#                    save_file=None, label_cut='', path_file=None):
#        
#        plt.figure(figsize=figsize)
#        results = {}
#        for j, key in enumerate(self.key_label):
#            if model_used==None:
#                results[key] = self.compute_residuals(filters, key_label=key, model_used=key, total_sn = True)
#                results[key] = results[key]['total']
#            else:
#                results[key] = self.compute_residuals(filters, key_label=key, model_used=model_used[j], total_sn = True)
#                results[key] = results[key]['total']
#
#        for i, f in enumerate(filters):
#            plt.subplot(len(filters), 1, i+1)
#            blab = f.replace(label_cut, '')
#            plt.ylabel('$\Delta$ '+blab)
#            for v, key in enumerate(self.key_label):
#                key_text = key.replace('_', ' ')
#                plt.hist(results[key][f]['residue'], bins=BINNING[i], normed=NORM, label=key_text, lw=3, histtype='step', color=color[v])
#                mean = np.mean(results[key][f]['residue'])
#                sig = np.std(results[key][f]['residue']) 
#                
#                plt.text(txt_pos_x[v], txt_pos_y[v], 'mean '+key_text+' = %.3f \n'%(mean)+' rms '+key_text+' = %.3f'%(sig))
#
#        plt.legend()
#        if path_file != None:
#            plt.savefig(path_file) 
#        plt.show()
#####################################################


    def res_by_phase(self, filters, phase, save_file=None, model_used=None,
                     color=['b', 'r', 'black', 'g'], figsize=(15, 25),
                     xlim=(-12,49), lim_res=(0.02,0.3),
                     label_cut='',path_file=None):
      
        mean = {}
        std = {}
        results = {}
        plt.figure(figsize=figsize)
        results = {}
        for s, key in enumerate(self.key_label):
            mean[key] = {}
            std[key] = {}
            if model_used==None:
                results[key] = self.compute_bin_residuals(filters, phase, key_label=key, model_used=key, total_sn = True)
                results[key] = results[key]['total']
            else:
                results[key] = self.compute_res(filters, phase, key_label=key, model_used=model_used[s], total_sn = True)
                results[key] = results[key]['total']
                
            for i, band in enumerate(filters):
                mean[key][band] = {'mean_value':[], 'err_mean':[]}
                std[key][band] = {'std_value':[], 'err_std':[], 'wrms_value':[],
                                   'wrms_err':[]}
                
        
        for i, band in enumerate(filters):
            p = []
            for j in range(len(phase)-1):
                    for key in self.key_label:     
                        res_phase = np.array(results[key][band]['res_bin_phase'][str(phase[j])+'_'+str(phase[j+1])])
                        var_phase = np.array(results[key][band]['var_bin_phase'][str(phase[j])+'_'+str(phase[j+1])])
                        err_mean = np.std(res_phase)/np.sqrt(len(res_phase)-1)
                        err_std = np.std(res_phase)/np.sqrt(2*len(res_phase))
                        mean[key][band]['mean_value'].append(np.mean(res_phase))
                        std[key][band]['std_value'].append(np.std(res_phase))
                        wrms = np.float(np.sqrt(np.sum(res_phase**2/var_phase)/np.sum(1/var_phase)))
#                        wrms = biweight_S(res_phase)
                        std[key][band]['wrms_value'].append(wrms)
                        std[key][band]['wrms_err'].append( np.sqrt(wrms*2*len(res_phase)/(2*np.sum(1/var_phase))))
                        mean[key][band]['err_mean'].append(err_mean)
                        std[key][band]['err_std'].append(err_std)
                    p.append((phase[j+1]+phase[j])/2)
    
            plt.subplot(len(filters)*2, 2, 2*i+1)
            for n, key in enumerate(self.key_label):
                key_text = key.replace('_', ' ')
                plt.errorbar(p, mean[key][band]['mean_value'], yerr=mean[key][band]['err_mean'], color=color[n], label=key_text)
                plt.plot([-15.,50.],[0.,0.], color='black')
                blab = band.replace(label_cut, '')
                plt.ylabel('$\Delta$ '+blab)
                plt.xlim(xlim)
                
            if i == 0:
                
                plt.title('Mean')
            if i == len(filters)-1:
                plt.xlabel('phase')
                plt.legend()    
            plt.subplot(len(filters)*2, 2, 2*(i+1))
#            plt.errorbar(p, std[key][band]['std_value'], yerr=std[key][band]['err_std'], color=color[n], label=key_text)
            for n, key in enumerate(self.key_label):
                key_text = key.replace('_', ' ')
#                plt.errorbar(p, std[key][band]['std_value'], yerr=std[key][band]['err_std'], color=color[n], label=key_text)
#            plt.subplot(len(filters), 1, i+1)
#            for n, key in enumerate(self.key_label):
#                key_text = key.replace('_', ' ')
                plt.errorbar(p, std[key][band]['wrms_value'], color=color[n], label=key_text)   
                plt.ylim(lim_res)
#                plt.errorbar(p, std[key][band]['std_value'], yerr=std[key][band]['err_std'], color=color[n], label=key_text)   
#                blab = band.replace(label_cut, '')
#                plt.ylabel('$\Delta$ '+blab)
                
                
            plt.xlim(xlim)
            if i == 0:
                plt.title('wRMS')
#                plt.legend()
            if i == len(filters)-1:
                plt.xlabel('phase')
                plt.legend()    
        if path_file != None:
            plt.savefig(path_file, bbox_inches='tight') 
        for key in self.key_label:
            pkl.dump(std[key], open(self.sad_path+self.sample+key+'wrms.pkl', 'wb'))
        plt.show()

    def plot_lc(self, sn_name, results, filters, scale_coeff=1.5, tlim= [-8, 48],
                color=['purple', 'b', 'y', 'r', 'k', 'g'],
                subplot_kargs = {'left':0.12, 'right':0.99, 'bottom':0.05, 'top':0.97, 'hspace':0.001},
                model_used='sugar',
                key_label='sugar',
                xlim=(-20,+55),
                ylim_gs0=(10,15),
                ylim_res=(-0.15,0.15),
                label_cut='', figsize=(15,20),
                units= 'mag',
                path_file=None):  
        
        #reload snls Filter 
        if self.sample == 'jla' :
                head, table = read_lc_jla(sn_name, sad_path=self.sad_path, model=model_used)      
            
        if len(color) <= len(filters):
            raise ValueError('color len should be equal or at least superior thant filter len')
            
        plt.figure(figsize=figsize)
        plt.subplots_adjust(**subplot_kargs)
        for band in filters:
            height_ratios=[]
            for i in range(len(filters)+1):
                if i ==0:
                    height_ratios.append(len(filters)+1)
                else:
                    height_ratios.append(1.5)
        dic_color = {}
       
        gs = gridspec.GridSpec(len(filters)+1, 1, height_ratios=height_ratios)
        for i, band in enumerate(filters):
            dic_color[band] = {'color': color[i],
                                'cst' : i*scale_coeff,
                                'gs' : gs[i+1]}


        self.sys = sncosmo.get_magsystem(self.sys_used)
        params = self.dic_res[key_label][sn_name]['res']['parameters']   
        model = self.build_model(params, model_used)
        time = np.linspace(tlim[0],tlim[1],50)
        for band in filters:
            
            #plot lc 
            plt.subplot(gs[0])
            sample = self.sample.replace('_', ' ')
            mod = model_used.replace('_', ' ')
            plt.title(sn_name + ' z = '+
                      str(self.dic_res[key_label][sn_name]['res']['parameters'][0])+
                      ' (' +mod+', '+sample+')' + ' $\chi^{2}$ = ' +
                      str(self.dic_res[key_label][sn_name]['res']['chisq']) + 
                      ' $\chi^{2}/ndof $ = ' +
                      str(self.dic_res[key_label][sn_name]['res']['chisq']/self.dic_res[key_label][sn_name]['res']['ndof']) )
            model_flux = model.bandflux(band, (time+params[1]),
                                        zp=results[sn_name][band]['zp'][0],
                                        zpsys=results[sn_name][band]['zpsys'][0] )
            
            if units == 'mag':
                model_mag = []
                for  fl in model_flux:
                    if fl >= 0.:
                        model_mag.append(self.sys.band_flux_to_mag(fl, band))
                    else:
                        model_mag.append(np.NaN)
                model_array = np.array(model_mag) + dic_color[band]['cst']
                obs = np.array(results[sn_name][band]['mag'])+dic_color[band]['cst']
                err_obs = np.array(results[sn_name][band]['data_err_mag'])
                
            else :            
                model_array = model_flux 
                obs = np.array(results[sn_name][band]['flux']) 
                err_obs = np.array(results[sn_name][band]['data_err_flux'])
            phase = np.array(results[sn_name][band]['phase'])
            
            plt.plot(time, model_array, 
                     dic_color[band]['color'],linestyle='--')
            for i, t in enumerate(phase):
                if t >= tlim[0] and  t <= tlim[1]:         
                    plt.errorbar(phase[i], obs[i],
                                 linestyle='', yerr=err_obs[i], xerr=None, mec=dic_color[band]['color'],
                                 color=dic_color[band]['color'], alpha=0.9, marker='o', zorder=0)
                else:
                    plt.errorbar(phase[i], obs[i],
                                 linestyle='', yerr=err_obs[i], xerr=None, mec=dic_color[band]['color'],
                                 color=dic_color[band]['color'], alpha=0.9, marker='x', zorder=0)                   
            plt.xlim(xlim) 
            plt.ylim(ylim_gs0)
            
            if units == 'mag':
                plt.ylabel('mag + const.')
                residuals = 'residue_mag'
            else:
                plt.ylabel('Flux')
                residuals = 'residue_flux'
                
            if units == 'mag' :
                plt.gca().invert_yaxis()               
            #plot residuals
            plt.subplot(dic_color[band]['gs'])
            residuals = np.array(results[sn_name][band][residuals])
            plt.errorbar(phase, residuals, yerr=err_obs, xerr=None, 
                         color=dic_color[band]['color'], 
                         alpha=0.9,
                         linestyle='', 
                         marker='o',
                         mec=dic_color[band]['color'],
                         zorder=0)
            plt.plot([-20,55], [0,0], dic_color[band]['color'], linestyle='--')
            plt.xlim(xlim)
            plt.ylim(ylim_res)
            for i, band in enumerate(filters):
                plt.subplot(gs[i+1])
                blab = band.replace(label_cut, '')
                plt.ylabel('$\Delta$ '+blab)
                if units == 'mag' :
                    plt.gca().invert_yaxis()
                if i == len(filters)-1 :
                    plt.xlabel('Phases (days)')

        if path_file != None:
            plt.savefig(path_file) 
        plt.show()

