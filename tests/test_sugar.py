#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  2 13:34:25 2018

@author: florian
"""
import cPickle
import numpy as np
import sncosmo
import read_data as RD
import builtins_SNF as Build_SNF
from matplotlib import pyplot as plt
from scipy import integrate
from sncosmo.salt2utils import BicubicInterpolator
import numpy as np

filters=['BSNf','VSNf','RSNf']
errorscale=True
width=20
###########
#constants
###########

CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27        # [erg s]

clight = 299792.458
H0 = 0.000070

try:
    Build_SNF.register_SNf_bands_width(width=width)
    Build_SNF.mag_sys_SNF_width(width=width)
    Build_SNF.register_SUGAR()    

except:
    print 'Filters and mag sys already registred' 



def test_plot_sugar():
    '''
    plot the difference betweem Mgr phot and Mgr spec
    '''
    
    SUGAR_parameter_pkl = '../sugar_model/sugar_parameters.pkl'
    dico = cPickle.load(open(SUGAR_parameter_pkl))
    Mgr_phot = []
    Mgr_spec = []
    x_zh = []
    list_sn = ['SN2006cj','SN2007kk']
    
#    for sn_name in list_sn:
    for sn_name in dico.keys():
        grey = dico[sn_name]['grey']
        q1 = dico[sn_name]['q1']
        q2 = dico[sn_name]['q2']
        q3 = dico[sn_name]['q3']
        av = dico[sn_name]['Av']
        cov_x = dico[sn_name]['cov_q']
        
        
        
        meta = cPickle.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
        data, zhl, mwebv, daymax = RD.read_meta_SNF(meta, sn_name, filters=filters, errorscale=errorscale, width=width)
        zcmb = meta[sn_name]['host.zcmb']
        
        def int_cosmo(z, Omega_M=0.3):     
            return 1./np.sqrt(Omega_M*(1+z)**3+(1.-Omega_M))
            
        def luminosity_distance(zhl,zcmb):
            
        
            integr = integrate.quad(int_cosmo, 0, zcmb)[0]
        
            return (1+zhl)*(clight/H0)*integr
         
        def distance_modulus_th(zhl,zcmb):      
            return 5.*np.log(luminosity_distance(zhl,zcmb))/np.log(10.)-5.                
        grey = grey + distance_modulus_th(zhl,zcmb)
        #num = 1
        #grey += num*5*np.log10(1+zhl)
        
        source = sncosmo.get_source('sugar')
        dust = sncosmo.CCM89Dust()
        model = sncosmo.Model(source=source, effects=[dust],effect_names=['mw'], effect_frames=['obs'])    
        model.set(mwebv=mwebv,z=zhl, q1=q1, q2=q2, q3=q3, A=av, Mgr=grey, t0=meta[sn_name]['salt2.DayMax'])
        
        #sncosmo.plot_lc(data, model=model)
        
        #res, fitted_model = sncosmo.fit_lc(data, model, ['t0'], modelcov = False, bounds={'t0':(daymax-5,daymax+5)},phase_range=(-12,42))
        
        res, fitted_model = sncosmo.fit_lc(data, model, ['Mgr'], modelcov = False)
        Mgr_spec.append(grey)
        Mgr_phot.append(res.parameters[6])
#        sncosmo.plot_lc(data, model=fitted_model,errors=res.errors)
        x_zh.append(5*np.log10(1+zhl))
#        print Mgr_spec, Mgr_phot
        #plt.savefig('../sugar_analysis_data/results/'+str(sn_name)+'_Mgr_fix_corr.pdf')
        plt.show()
    
    Mgr_spec = np.array(Mgr_spec)
    Mgr_phot = np.array(Mgr_phot)
    diff_Mgr = Mgr_phot - Mgr_spec
    plt.plot(x_zh,diff_Mgr,'o')
    plt.xlabel('5*log_10(1+zhl)')
    plt.ylabel('Mgr_phot - Mgr_spec')
    plt.savefig('../sugar_analysis_data/results/dep_Mgr_phot_spec.pdf')
    print np.corrcoef(x_zh,diff_Mgr)[0,1]
    plt.show()
    
def plot_sugar_model():
    """
    plot the point of sugar model and this interpotion 
    to this if the interpotion is good
    """
    source = sncosmo.get_source('sugar')
    dust = sncosmo.CCM89Dust()
    model = sncosmo.Model(source=source)  
    sugarm = np.genfromtxt('../sugar_model/SUGAR_model_v1.asci')
#    print sugarm
    M0_5613 = []
    M1_5613 = []
    M2_5613 = []
    M3_5613 = []
    M4_5613 = []
    phase_5613 = []
    for i in range(len(sugarm)):
    #    print type(wave[i])
        if sugarm[i,1] == 5632.30974:
            M0_5613.append(sugarm[i,2])
            M1_5613.append(sugarm[i,3])
            M2_5613.append(sugarm[i,4])
            M3_5613.append(sugarm[i,5])
            M4_5613.append(sugarm[i,6])
            phase_5613.append(sugarm[i,0])
    M0_5613 = np.array(M0_5613)
    M1_5613 = np.array(M1_5613)
    M2_5613 = np.array(M2_5613)
    M3_5613 = np.array(M3_5613)
    M4_5613 = np.array(M4_5613)
#    flux = (10. ** (-0.4 * (M0_5613 +   M1_5613+  M2_5613 +  M3_5613 +  M4_5613 + 40. + 48.59)) / (5632.30974 ** 2 / 299792458. * 1.e-10))
    #print   M0_5613     
#    plt.plot(phase_5613,flux,'ro')
    time=np.linspace(-20,63,100)
    plt.plot(time,model.flux(np.linspace(-20,63,100),np.array([7000])))
#    plt.plot(time,model.flux(np.linspace(-20,50,100),np.array([5632.30974])))
    plt.xlabel('Phase')
    plt.ylabel('Flux')
    plt.text(20,6e-18,'Wave = 5632.30974 A')
    plt.legend()
    pdffile = '../sugar_analysis_data/results/sugar_model.pdf'
    plt.savefig(pdffile, bbox_inches='tight')  
    plt.show()    