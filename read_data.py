# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 17:50:29 2017

@author: mondon
"""
import cPickle
import numpy as np
#import sugar
from scipy import integrate
from astropy.table import Table
import sncosmo
from math import isnan
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
import pyfits
import sys
from math import log, log10
###########
#constants
###########

CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27        # [erg s]

clight = 299792.458
H0 = 0.000070
###############
# wavelength limits for salt2 model
wl_min_sal = 3000
wl_max_sal = 7000

# wavelength limits for sugar model
wl_min_sug = 3341.41521
wl_max_sug = 8576.61898

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
    print len(salt_parm)
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

def wl_cut_salt2(fname, EBV, z):
    # wavelength limits for salt2 model
    wl_min_sal = 3000
    wl_max_sal = 7000
    

    filt = sncosmo.get_bandpass(fname)
    wlen = filt.wave
    tran = filt.trans
    dt = 10000
    spl = Spline1d(wlen, tran, k=1, ext = 1)

    xs = np.linspace(min(wlen), max(wlen), dt)
    dxs = ((max(wlen)-min(wlen))/(dt-1))
    wlen_eff = np.sum((10**(-A_l(xs,EBV)/2.5))*spl(xs)*xs*dxs)/np.sum((10**(-A_l(xs,EBV)/2.5))*spl(xs)*dxs)
    if wl_min_sal >= wlen_eff/(1+z) or wlen_eff/(1+z) >= wl_max_sal:
        return ('False', wlen_eff/(1+z))
    else:
        return('True', wlen_eff/(1+z))

def wl_cut_sugar(fname, z):		
        filt = sncosmo.get_bandpass(fname)
        wlen = filt.wave
        tran = filt.trans
        dt = 10000
        wlen_shift = wlen/(1+z)
        spl = Spline1d(wlen_shift, tran, k=1, ext = 1)
        xs = np.linspace(min(wlen_shift), max(wlen_shift), dt)
        dxs = ((max(wlen_shift)-min(wlen_shift))/(dt-1)) 
        area_full = np.sum(spl(xs)*dxs) # full area under the filter

        xs_cut = np.linspace(wl_min_sug, wl_max_sug, dt)
        dxs_cut = ((wl_max_sug-wl_min_sug)/(dt-1))
        area_cut = np.sum(spl(xs_cut)*dxs_cut)
        r = 1.-area_cut/area_full # area under the filter outside the model

        if r < 0.1:
            return ('True',r)
        else:
            return ('False',r)





def A_l(xx, EBV, Rv=3.1):
    A_l = []
    xx = np.array(xx)
    x = 1./(xx/10**4.) # xx in AA
    for i in x:
        if i < 1.1:
            y = i**1.61
            a = 0.574*y
            b = -0.527*y
        elif 1.1 <= i <= 3.3:
            y = i - 1.82
            a = 1 + y*(0.17699 + y*(-0.50447 + y*(-0.02427+ y*(0.72085 + y*(0.01979 + y*(-0.77530 + 0.32999*y))))))
            b =  y*(1.41338 + y*(2.28305 + y*(1.07233 +y*(-5.38434 + y*(-0.62251 + y*(5.30260 - 2.09002*y))))))
        elif 3.3 < i < 5.9:
            a = 1.752 - 0.316*i - 0.104/((i-4.67)*(i-4.67) + 0.341)
            b = -3.09 + 1.825*i + 1.206/((i-4.62)*(i-4.62) + 0.263)
        elif 5.9 <= i <= 8.:
            Fa = -(i-5.9)*(i-5.9) * (0.04473 + 0.009779*(i-5.9))
            a = 1.752 - 0.316*i - 0.104/((i-4.67)*(i-4.67) + 0.341) + Fa
            Fb = (i-5.9)*(i-5.9) * (0.213 + 0.1207*(i-5.9))
            b = -3.09 + 1.825*i + 1.206/((i-4.62)*(i-4.62) + 0.263) + Fb
        A = (a*Rv + b)*EBV
        A_l.append(A)
    A_l = np.array(A_l)
    return A_l


            
def mag_to_flux(mag,band,width):
    vega = sncosmo.get_magsystem('vega_snf_'+str(width))
    flux = vega.band_mag_to_flux(mag, band)

    return flux
   
def get_zp(band):
        if band == 'USNf' :
            band_file = 'USNf_3300-4102.dat'
            
        if band == 'BSNf' :
            band_file = 'BSNf_4102-5100.dat'
            
        if band == 'VSNf' :
            band_file = 'VSNf_5200-6289.dat'
            
        if band == 'RSNf' :
            band_file = 'RSNf_6289-7607.dat'
            
        if band == 'ISNf' :
            band_file = 'ISNf_7607-9200.dat'
            
            
        filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file)
        wlen = filt2[:,0]
        tran = filt2[:,1]
        splB = Spline1d(wlen, tran, k=1,ext = 1)    
#        
        fits = pyfits.open('../sugar_analysis_data/Vega.fits')
        fit = fits[1]
        wl = np.zeros(len(fit.data))
        flux = np.zeros(len(fit.data))
        for i in range(len(fit.data)):
             wl[i] = fit.data[i][0]
             flux[i] = fit.data[i][1]
        splref = Spline1d(wl, flux, k=1,ext = 1)         
    
#        band_sn = sncosmo.get_bandpass(band)
#        wlen = band_sn.wave
#        tran = band_sn.trans
#        splB = Spline1d(wlen, tran, k=1,ext = 1)    
        
        #computation of the integral
        dt = 100000
        xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
        dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))
        
        
        I2 = np.sum(splref(xs)*xs/ (HPLANCK * CLIGHT)*splB(xs)*dxs)       

        return 2.5*log10(I2)    
        
def read_meta_SNF(meta,sn_name,filters=['BSNf','VSNf','RSNf'],model='salt2',errorscale=True, width=20):
    """
    """
    snfit_test = open('snfit_test.txt','w')
    # wavelength limits for salt2 model
    wl_min_sal = 3000
    wl_max_sal = 7000
    
    # wavelength limits for sugar model
    wl_min_sug = 3341.41521
    wl_max_sug = 8576.61898
    vega = sncosmo.get_magsystem('vega_snf_'+str(width))
    time = []
    band = []
    flux = []
    fluxerr = []
    zp = []
    zpsys = []
    h = 10**-9
    errorscale_factor = meta[sn_name]['target.errorscale']
    zhl = meta[sn_name]['host.zhelio']
    zcmb = meta[sn_name]['host.zcmb']
    daymax = meta[sn_name]['salt2.DayMax']
    mwebv = meta[sn_name]['target.mwebv']
    sn_data = meta[sn_name]['spectra']
    for t in sn_data.keys():

        for f in filters:
            if f == 'USNf'or f == 'BSNf':
                try:
                    quality_flag = sn_data[t]['procB.Quality']
                    
                except:
                    quality_flag = 1
                    
            elif f == 'VSNf'or f == 'RSNf' or f == 'ISNf':
                try:
                    quality_flag = sn_data[t]['procR.Quality']
                except:
                    quality_flag = 1
            else:
                print >> sys.stderr, \
                    'filter not recognise' 
                return
                
#            if isnan(sn_data[t]['mag.'+f]) == False:
            if quality_flag == 1 and not isnan(sn_data[t]['mag.'+f]):
#
                time.append(sn_data[t]['obs.mjd'])
                
                if f == 'USNf':
                    fn = 'fU_'+str(width)
                    band.append(fn)
                if f == 'BSNf':
                    fn = 'fB_'+str(width)
                    band.append(fn)
                if f == 'VSNf':
                    fn = 'fV_'+str(width)
                    band.append(fn)
                if f == 'RSNf':
                    fn = 'fR_'+str(width)
                    band.append(fn)
                if f == 'ISNf':
                    fn = 'fI_'+str(width)
                    band.append(fn) 
#                if f == 'BSNf':
#                    band.append('fB')
#                if f == 'VSNf':
#                    band.append('fV')
#                if f == 'RSNf':
#                    band.append('fR')                   
#                band.append(f)
                
                flux.append(mag_to_flux(sn_data[t]['mag.'+f],fn,width))
                zp.append(2.5*np.log10(vega.zpbandflux(fn)))


#                print '~~~~~~~~~~~~~~~~~~~'
#                print 2.5*np.log10(vega.zpbandflux(f))
#                print get_zp(f)
#                print '~~~~~~~~~~~~~~~~~~~'
#                zp.append(get_zp(f))
                
                zpsys.append('vega_snf_'+str(width))
                if errorscale:

                    err_mag = sn_data[t]['mag.'+f+'.err']

                else:  
                    err_mag = sn_data[t]['mag.'+f+'.err']/errorscale_factor

#                fluxerr.append(((np.abs(mag_to_flux(sn_data[t]['mag.'+f]+h,f)-mag_to_flux(sn_data[t]['mag.'+f] ,f)) / h))*  err_mag)
                fluxerr.append(err_mag * mag_to_flux(sn_data[t]['mag.'+f],fn,width) / 1.0857362047581294)
#                snfit_test.write(str(sn_data[t]['obs.mjd'])+' '+str(f)+' '+str(mag_to_flux(sn_data[t]['mag.'+f],f))+' '+str(err_mag * mag_to_flux(sn_data[t]['mag.'+f],f) / 1.0857362047581294)+' '+str(2.5*np.log10(vega.zpbandflux(f)))+' '+'Vega'+'\n')
    data = Table([time, band, flux, fluxerr, zp, zpsys], names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), meta={'name': 'data'})


    
    snfit_test.close()
    dic = {}
    for x in data:
        if x[1] in dic.keys():
            dic[x[1]] += 1
        else:
            dic[x[1]] = 1
            
    f_in = {}	
    f_out = {}
    for i in dic.keys():
        if model == 'salt2':
            res = wl_cut_salt2(i, mwebv, zhl) 
            if res[0] == 'True':
                f_in[i] = res[1]
            else:
                f_out[i] = res[1]
                print 'We excluded passband %s (%d points) because restframewavelength = %7.3f does not belong to the interval [%d,%d]' % (i,dic[i],res[1],wl_min_sal,wl_max_sal)
        elif model == 'sugar':
            res = wl_cut_sugar(i, zhl)
            if res[0] == 'True':
                f_in[i] = res[1]
            else:
                f_out[i] = res[1]
                print 'We excluded passband %s (%d points) because it does not belong to the interval [%d,%d]' % (i,dic[i],wl_min_sug,wl_max_sug)
        else:
            print 'ERROR: model name has to be salt2 or sugar'

    mask = []
    for row in data:
        if row[1] in f_in.keys():
            mask.append(True)
        else:
            mask.append(False)
    mask = np.array(mask)

    data_cut = sncosmo.select_data(data, mask)
    return data_cut, zhl, mwebv, daymax
    
    
def mag_to_flux_hand(mag,band):
    #interpolation of TB and Trest
    if band == 'USNf' :
        band_file = 'USNf_3300-4102.dat'

    if band == 'BSNf' :
        band_file = 'BSNf_4102-5100.dat'

    if band == 'VSNf' :
        band_file = 'VSNf_5200-6289.dat'

    if band == 'RSNf' :
        band_file = 'RSNf_6289-7607.dat'
    
    if band == 'ISNf' :
        band_file = 'ISNf_7607-9200.dat'
       
        
    filt2 = np.genfromtxt('../sugar_analysis_data/data/Instruments/SNIFS/'+ band_file)
    wlen = filt2[:,0]
    tran = filt2[:,1]
    splB = Spline1d(wlen, tran, k=1,ext = 1)

#    fits = pyfits.open('../sugar_analysis_data/Vega.fits')
#    fit = fits[1]
#    wl = np.zeros(len(fit.data))
#    flux = np.zeros(len(fit.data))
#    for i in range(len(fit.data)):
#         wl[i] = fit.data[i][0]
#         flux[i] = fit.data[i][1]
#    splref = Spline1d(wl, flux, k=1,ext = 1)

    #interpolation of ref spectrum
    data = np.genfromtxt('../sugar_analysis_data/data/MagSys/bd_17d4708_stisnic_002.ascii')
    dispersion = data[:,0]
    flux_density = data[:,1]  
    splref = Spline1d(dispersion, flux_density, k=1,ext = 1)
#    
    
    
    #computation of the integral
    dt = 100000
    xs = np.linspace(float(wl_min_sal), float(wl_max_sal), dt)
    dxs = (float(wl_max_sal-wl_min_sal)/(dt-1))
    
#    vega = sncosmo.get_magsystem('vega_snf')
#    print vega.band_mag_to_flux(mag, band)
#    print '##################################################'
    I2 = np.sum(splref(xs)*xs/ (HPLANCK * CLIGHT)*splB(xs)*dxs)
#    print vega.zpbandflux(band)
#    print I2
#    print '##################################################'
    return I2*10**(-0.4*(mag))
    
    
    
    
    
    
    
    
    
    