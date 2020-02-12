#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 16:09:44 2019

@author: florian
"""
import pylab as plt
import numpy as np
import sncosmo
import matplotlib.gridspec as gridspec
try:
   import cPickle as pkl
except ModuleNotFoundError:
   import pickle as pkl


from .builtins import register_SNf_bands_width, mag_sys_SNF_width,  builtins_jla_bandpasses, mag_sys_jla

def plot_lc_res(sn_name, 
            dic_path='../../sugar_analysis_data/resfitlc_csp_sugar.pkl', 
            sed_model = 'sugar',
            sample='csp',
            filters = ['cspb', 'cspg', 'cspv3014', 'cspv3009', 'cspv9844', 'cspr', 'cspi'],
            outfile_path='../../sugar_analysis_data/results/save_plot/jla/'):
    """
    doc a faire 
    """
    if sample == 'jla':
        builtins_jla_bandpasses()
        mag_sys_jla()
        sample='jla_Vega2'
    elif sample == 'SNf':
        register_SNf_bands_width(width=10)
        mag_sys_SNF_width(width=10)
        sample = 'vega_snf_10'
        
    dic = pkl.load(open(dic_path))
    sn_title = ''
    filt =[]
    if type(dic['data'][sn_name]['res']) is str:
        return 'Fit fail for '+sn_name
    for f in filters:
#        minwave = (sncosmo.get_bandpass(f)).minwave()
        if f in dic['data'][sn_name]['data_table']['band']: #and minwave/(1+float(dic['data'][sn_name]['res']['parameters'][0]))>=3254.01639:
            filt.append(f)
    for letter in sn_name:
        if letter != '_' or letter !='.':
            sn_title += letter
        else:
            break
    
    sys = sncosmo.get_magsystem(sample)
    
    if sed_model is 'sugar':
        source = sncosmo.get_source('sugar')
        dust = sncosmo.CCM89Dust()
        model = sncosmo.Model(source=source,
                              effects=[dust],
                              effect_names=['mw'],
                              effect_frames=['obs'])

        model.set(z = float(dic['data'][sn_name]['res']['parameters'][0]))
        model.set(t0 = float(dic['data'][sn_name]['res']['parameters'][1]))
        model.set(Xgr = float(dic['data'][sn_name]['res']['parameters'][2]))
        model.set(q1 = float(dic['data'][sn_name]['res']['parameters'][3]))
        model.set(q2 = float(dic['data'][sn_name]['res']['parameters'][4]))
        model.set(q3 = float(dic['data'][sn_name]['res']['parameters'][5]))
        model.set(A = float(dic['data'][sn_name]['res']['parameters'][6]))
        model.set(mwebv = float(dic['data'][sn_name]['res']['parameters'][7]))
        
        
    if sed_model is 'salt2':
        source = sncosmo.get_source('salt2')
        dust = sncosmo.CCM89Dust()
        model = sncosmo.Model(source=source,
                              effects=[dust],
                              effect_names=['mw'],
                              effect_frames=['obs'])

        model.set(z = float(dic['data'][sn_name]['res']['parameters'][0]))
        model.set(t0 = float(dic['data'][sn_name]['res']['parameters'][1]))
        model.set(x0 = float(dic['data'][sn_name]['res']['parameters'][2]))
        model.set(x1 = float(dic['data'][sn_name]['res']['parameters'][3]))
        model.set(c = float(dic['data'][sn_name]['res']['parameters'][4]))
        model.set(mwebv = float(dic['data'][sn_name]['res']['parameters'][5]))
    
    model_title = sed_model
    time = np.linspace(-12,48,50)
    t0 = float(dic['data'][sn_name]['res']['parameters'][1])
    
    for f in filt :
        Filtre = dic['data'][sn_name]['data_table']['band'] == f
        Filtre &= dic['data'][sn_name]['data_table']['time'] > -12 + t0
        Filtre &= dic['data'][sn_name]['data_table']['time'] < 48 + t0

        flux = dic['data'][sn_name]['data_table']['flux'][Filtre]
        flux_err = dic['data'][sn_name]['data_table']['fluxerr'][Filtre]
        time_obs = dic['data'][sn_name]['data_table']['time'][Filtre] - t0
        band_obs = dic['data'][sn_name]['data_table']['band'][Filtre]
        zp = float(dic['data'][sn_name]['data_table']['zp'][Filtre][0])
        try: 
            mod_flux = model.bandflux(f, time+t0, zp=zp, zpsys=sample)
        except:  
            filt.remove(f) 
            
    plt.figure(figsize=(7,12))
    plt.subplots_adjust(left=0.12, right=0.99,bottom=0.05,top=0.97,hspace=0.001)
    
    
    height_ratios=[]
    for i in range(len(filt)+1):
        if i ==0:
            height_ratios.append(len(filt)+1)
        else:
            height_ratios.append(1.5)

    gs = gridspec.GridSpec(len(filt)+1, 1,height_ratios=height_ratios)
    dic_color = {}
    col = ['purple', 'b', 'y', 'r', 'k', 'g']
    
    for i, band in enumerate(filt):
        dic_color[band] = {'color': col[i],
                            'cst' : i,
                            'gs' : gs[i+1]}
        
    

    
    dic_residual = {}
    for i, band in enumerate(filt):
        dic_residual[band] = None
    

            
    for f in filt:
        

        #print f
        Filtre = dic['data'][sn_name]['data_table']['band'] == f
        Filtre &= dic['data'][sn_name]['data_table']['time'] > -12 + t0
        Filtre &= dic['data'][sn_name]['data_table']['time'] < 48 + t0

        flux = dic['data'][sn_name]['data_table']['flux'][Filtre]
        flux_err = dic['data'][sn_name]['data_table']['fluxerr'][Filtre]
        time_obs = dic['data'][sn_name]['data_table']['time'][Filtre] - t0
        band_obs = dic['data'][sn_name]['data_table']['band'][Filtre]
        zp = float(dic['data'][sn_name]['data_table']['zp'][Filtre][0])

                
        mag_err_obs = []
        mag_obs = []
        for i, fl in enumerate(flux):
            if fl >= 0.:
                mag_obs.append(sys.band_flux_to_mag(fl, f))
                mag_err_obs.append((flux_err[i]/fl) * 1.0857362047581294)
            else:
                mag_obs.append(np.NaN)
                mag_err_obs.append(np.NaN)
                
        mag_err_obs = np.array(mag_err_obs)
        mag_obs = np.array(mag_obs)

      
        mod_flux = model.bandflux(f, time+t0, zp=zp, zpsys=sample)   #zp is the same for the same filter
        mod_flux_obs = model.bandflux(band_obs, time_obs+t0, zp=zp, zpsys=sample) 
        
        mag = []
        for  fl in mod_flux:
            if fl >= 0.:
                mag.append(sys.band_flux_to_mag(fl, f))
            else:
                mag.append(np.NaN)
        mag = np.array(mag)
        
        mag_phase_obs = []
        for  fl in mod_flux_obs:
            if fl >= 0.:
                mag_phase_obs.append(sys.band_flux_to_mag(fl, f))
            else:
                mag_phase_obs.append(np.NaN)
        mag_phase_obs = np.array(mag_phase_obs)
        




        plt.subplot(gs[0])
        plt.plot(time, mag+dic_color[f]['cst'], dic_color[f]['color'],linestyle='--')
        plt.scatter(time_obs, mag_obs+dic_color[f]['cst'], c=dic_color[f]['color'])
        plt.errorbar(time_obs, mag_obs+dic_color[f]['cst'],
                     linestyle='', yerr=mag_err_obs, xerr=None, 
                     ecolor=dic_color[f]['color'],alpha=0.9,marker='.',zorder=0)
        
    
        if f == filt[len(filt)-1]: 
            YMAX = np.max(mag_obs+dic_color[f]['cst']) + 1.
            if type(YMAX) != float: 
                YMAX = 22.
        if f == filt[0]:
            YMIN = np.min(mag_obs+dic_color[f]['cst']) - 1.
            if type(YMIN) != float: 
                YMIN = 15.
            
        plt.subplot(dic_color[f]['gs'])
        plt.scatter(time_obs, mag_obs - mag_phase_obs, c=dic_color[f]['color'])
        plt.plot([-15,55], [0,0], dic_color[f]['color'], linestyle='--')
        plt.xlim(-15,+55)
        plt.ylim(-0.15,0.15)
        
        dic_residual[f] = mag_obs - mag_phase_obs

    plt.subplot(gs[0])
    plt.xlim(-15,+55)
    plt.ylim(YMIN,YMAX)
    plt.title(sn_title + ' (' +model_title+', '+sample+')', fontsize=18)
    plt.subplot(gs[0])
    plt.ylabel('mag + const.', fontsize=18)
    plt.gca().invert_yaxis()
    if sample == 'vega_snf_10' :
        filt = ['USNf', 'BSNf', 'VSNf', 'RSNf', 'ISNf']
    for i, band in enumerate(filt):
        plt.subplot(gs[i+1])
        if '_' in band:
            band = band.split('_')[1]
            band.replace('::', ' ')
        plt.ylabel('$\Delta$ '+band , fontsize=18)
        plt.gca().invert_yaxis()
        if i == len(filt)-1 :
            plt.xlabel('Phases (days)', fontsize=18)
            plt.xlim(-15,+55)
    
    plt.savefig(outfile_path+'%s_%s.png'%((sn_title,model_title)))
#    plt.show()
    plt.close()
    return dic_residual