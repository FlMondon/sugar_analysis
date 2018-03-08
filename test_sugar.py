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




SUGAR_parameter_pkl = '../sugar/sugar/data_output/data_output/sugar_parameters.pkl'
dico = cPickle.load(open(SUGAR_parameter_pkl))
sn_name = 'SN2007kk'
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


source = sncosmo.get_source('sugar')
dust = sncosmo.CCM89Dust()
model = sncosmo.Model(source=source, effects=[dust],effect_names=['mw'], effect_frames=['obs'])    
model.set(mwebv=mwebv,z=zhl, q1=q1, q2=q2, q3=q3, A=av, Mgr=grey, t0=meta[sn_name]['salt2.DayMax'])
res, fitted_model = sncosmo.fit_lc(data, model, ['t0'], modelcov = False, bounds={'t0':(daymax-5,daymax+5)},phase_range=(-12,42))
sncosmo.plot_lc(data, model=fitted_model)
plt.savefig('../sugar_analysis_data/results/'+str(sn_name)+'_t0_free_between5days.pdf')
plt.show()

sugarm = np.genfromtxt('../sugar_analysis_data/SUGAR_model_v1.asci')
M0_5613 = []
M1_5613 = []
M2_5613 = []
M3_5613 = []
M4_5613 = []
phase_5613 = []
for i in range(len(sugarm)):
#    print type(wave[i])
    if sugarm[i,1] == 5613.07529000:
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
flux = (10. ** (-0.4 * (M0_5613 +   M1_5613+  M2_5613 +  M3_5613 +  M4_5613 + 40. + 48.59)) / (5613.07529000 ** 2 / 299792458. * 1.e-10))
#print   M0_5613     
#plt.plot(phase_5613,flux,'ro')
#time=np.linspace(-12,42,100)
#plt.plot(time,model.flux(np.linspace(-12,42,100),np.array([5613.07529000])))