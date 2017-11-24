# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 17:50:29 2017

@author: mondon
"""
import cPickle
import numpy as np
import sugar
from scipy import integrate
from astropy.table import Table
import sncosmo

###########
#constants
###########

clight = 299792.458
H0 = 0.000070
###############


def read_sugar_salt2_parameters():
    """
    """
    SUGAR_parameter_pkl = '../sugar/sugar/data_output/data_output/sugar_parameters.pkl'
    lds = sugar.load_data_sugar()
    lds.load_salt2_data()
    Filtre = np.array([True]*len(lds.X0))
    mb = lds.mb
    mb_err = lds.mb_err
    x1 = lds.X1
    x1_err = lds.X1_err
    c = lds.C
    c_err = lds.C_err
    cov_mb_x1 = lds.X1_mb_cov
    cov_mb_C = lds.C_mb_cov
#    cov_mb_x1 = np.zeros(len(lds.C_mb_cov))
#    cov_mb_C = np.zeros(len(lds.C_mb_cov))
    cov_x1_C = lds.X1_C_cov  
#    print cov_mb_x1
    grey = np.zeros_like(mb)
    q1 = np.zeros_like(mb)
    q2 = np.zeros_like(mb)
    q3 = np.zeros_like(mb)
    av = np.zeros_like(mb)
    cov_x = np.zeros((len(mb),5,5))
    zcmb = lds.zcmb
    zhl = lds.zhelio
    zerr = lds.zerr
    dico = cPickle.load(open(SUGAR_parameter_pkl))
    
    
    for i in range(len(lds.sn_name)):
        sn = lds.sn_name[i]
        if sn in dico.keys():
             grey[i] = dico[sn]['grey']
             q1[i] = dico[sn]['q1']
             q2[i] = dico[sn]['q2']
             q3[i] = dico[sn]['q3']
             av[i] = dico[sn]['Av']
             cov_x[i] = dico[sn]['cov_q']
        else:
            Filtre[i] = False
    
    mb = mb[Filtre]
    mb_err = mb_err[Filtre]

    x1_err = x1_err[Filtre]
    c_err = c_err[Filtre]
    x1 = x1[Filtre]
    c = c[Filtre]

    cov_mb_x1 = cov_mb_x1[Filtre]
    cov_mb_C = cov_mb_C[Filtre]
    cov_x1_C = cov_x1_C[Filtre]
    cov_y = np.zeros((len(mb)*3,len(mb)*3))
    
    for i in range (len(mb)):
        cov_y[i*3,i*3] = mb_err[i]**2
        cov_y[i*3+ 1,i*3+ 1] = x1_err[i]**2
        
        cov_y[i*3+ 2,i*3+ 2] = c_err[i]**2
        cov_y[i*3+ 0,i*3+ 1] = cov_mb_x1[i]
        cov_y[i*3+ 1,i*3+ 0] = cov_mb_x1[i]
        cov_y[i*3+ 0,i*3+ 2] = cov_mb_C[i]
        cov_y[i*3+ 2,i*3+ 0] = cov_mb_C[i]
        cov_y[i*3+ 1,i*3+ 2] = cov_x1_C[i]      
        cov_y[i*3+ 2,i*3+ 1] = cov_x1_C[i]
    cov_salt = cov_y
    
    zcmb = zcmb[Filtre]
    zhl = zhl[Filtre]
    zerr = zerr[Filtre]
    cov_x = cov_x[Filtre]
    grey = grey[Filtre]
    q1 = q1[Filtre]
    q2 = q2[Filtre]
    q3 = q3[Filtre]
    av = av[Filtre]
    
    def int_cosmo(z, Omega_M=0.3):     
        return 1./np.sqrt(Omega_M*(1+z)**3+(1.-Omega_M))
        
    def luminosity_distance(zhl,zcmb):
        

        integr = np.zeros_like(zcmb)
        for i in range(len(zcmb)):
            integr[i] = integrate.quad(int_cosmo, 0, zcmb[i])[0]
    
        return (1+zhl)*(clight/H0)*integr
 
    def distance_modulus_th(zhl,zcmb):      
        return 5.*np.log(luminosity_distance(zhl,zcmb))/np.log(10.)-5.                
    grey = grey + distance_modulus_th(zhl,zcmb)               
    sug_parm = np.array([grey,q1,q2,q3,av]).T
    salt_parm = np.array([mb,x1,c]).T    
    
    cov_mat = np.zeros([len(grey)*5, len(grey)*5])
    cv = cov_x
    for i in range (len(grey)):
        cov_mat[i*5,i*5] = cv[i,0,0]
        cov_mat[i*5 +1,i*5] = cv[i,1,0]
        cov_mat[i*5 +2,i*5] = cv[i,2,0]
        cov_mat[i*5 +3,i*5] = cv[i,3,0]
        cov_mat[i*5 +4,i*5] = cv[i,4,0]
        cov_mat[i*5,i*5 +1] = cv[i,0,1]
        cov_mat[i*5 +1,i*5 +1] = cv[i,1,1]
        cov_mat[i*5 +2,i*5 +1] = cv[i,2,1]
        cov_mat[i*5 +3,i*5 +1] = cv[i,3,1]
        cov_mat[i*5 +4,i*5 +1] = cv[i,4,1]
        cov_mat[i*5 ,i*5 +2] = cv[i,0,2]
        cov_mat[i*5 +1,i*5 +2] = cv[i,1,2]
        cov_mat[i*5 +2,i*5 +2] = cv[i,2,2]
        cov_mat[i*5 +3,i*5 +2] = cv[i,3,2]
        cov_mat[i*5 +4,i*5 +2] = cv[i,4,2]
        cov_mat[i*5,i*5 +3] = cv[i,0,3]
        cov_mat[i*5 +1,i*5 +3] = cv[i,1,3]
        cov_mat[i*5 +2,i*5 +3] = cv[i,2,3]
        cov_mat[i*5 +3,i*5 +3] = cv[i,3,3]
        cov_mat[i*5 +4,i*5 +3] = cv[i,4,3]
        cov_mat[i*5,i*5 +4] = cv[i,0,4]
        cov_mat[i*5 +1,i*5 +4] = cv[i,1,4]
        cov_mat[i*5 +2,i*5 +4] = cv[i,2,4]
        cov_mat[i*5 +3,i*5 +4] = cv[i,3,4]
        cov_mat[i*5 +4,i*5 +4] = cv[i,4,4]
    
    cov_mat_grey = np.zeros([len(grey), len(grey)])

    for i in range (len(grey)):
        cov_mat_grey[i,i] = cv[i,0,0]
       


    
    return sug_parm, salt_parm, cov_salt, cov_mat_grey, cov_mat,  zhl, zcmb, zerr
    
def read_UBVRI():
    """
    """
    pkl_file = '../sugar_analysis_data/UBVRI_unix.pkl'
    dico = cPickle.load(open(pkl_file)) 
    zhl = []
    zcmb = []
    zerr = []
    mb = []
    x1 = []
    c = []
    mb_err = []
    x1_err = []
    c_err = []
    cov_mb_C = []
    cov_mb_x1 = []
    cov_x1_C = []
    sn_type = []
    for sn_name in dico.keys():
        zhl.append(dico[sn_name]['host.zhelio'])
        zcmb.append(dico[sn_name]['host.zcmb'])
        zerr.append(dico[sn_name]['host.zhelio.err'])
        mb.append(dico[sn_name]['salt2.RestFrameMag_0_B'])
        x1.append(dico[sn_name]['salt2.X1'])
        c.append(dico[sn_name]['salt2.Color']) 
        mb_err.append(dico[sn_name]['salt2.RestFrameMag_0_B.err'])
        x1_err.append(dico[sn_name]['salt2.X1.err'])
        c_err.append(dico[sn_name]['salt2.Color.err'])
        cov_mb_C.append(dico[sn_name]['salt2.CovColorRestFrameMag_0_B'])
        cov_mb_x1.append(dico[sn_name]['salt2.CovRestFrameMag_0_BX1'])
        cov_x1_C.append(dico[sn_name]['salt2.CovColorX1'])
        sn_type.append(dico[sn_name]['target.type'])
    
    cov_y = np.zeros((len(mb)*3,len(mb)*3))
    
    for i in range (len(mb)):
        cov_y[i*3,i*3] = mb_err[i]**2
        cov_y[i*3+ 1,i*3+ 1] = x1_err[i]**2
        
        cov_y[i*3+ 2,i*3+ 2] = c_err[i]**2
        cov_y[i*3+ 0,i*3+ 1] = cov_mb_x1[i]
        cov_y[i*3+ 1,i*3+ 0] = cov_mb_x1[i]
        cov_y[i*3+ 0,i*3+ 2] = cov_mb_C[i]
        cov_y[i*3+ 2,i*3+ 0] = cov_mb_C[i]
        cov_y[i*3+ 1,i*3+ 2] = cov_x1_C[i]      
        cov_y[i*3+ 2,i*3+ 1] = cov_x1_C[i]
    cov_salt = cov_y
    salt_parm = np.array([mb,x1,c]).T
    zhl = np.array(zhl)
    zcmb = np.array(zcmb)
    zerr = np.array(zerr)
    
    return salt_parm ,cov_salt, zhl, zcmb, zerr
    
def mag_to_flux(mag,band):
    vega = sncosmo.get_magsystem('vega')
    flux = vega.band_mag_to_flux(mag, band)    
    return flux
    
def read_meta_SNF(meta,sn_name,filters=['BSNf','VSNf','RSNf']):
    """
    """
    time = []
    band = []
    flux = []
    fluxerr = []
    zp = []
    zpsys = []
    sn_data = meta[sn_name]['spectra']
    for t in sn_data.keys():
        for f in filters:
            time.append(sn_data[t]['obs.mjd'])
            band.append(f)
            flux.append(mag_to_flux(sn_data[t]['mag.'+f],f))
            flux.append(mag_to_flux(sn_data[t]['mag.'+f+'err'],f))
            zp.append(0.)
            zpsys.append('vega')
            
    data = Table([time, band, flux, fluxerr, zp, zpsys], names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), meta={'name': 'data'})
    return data
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    