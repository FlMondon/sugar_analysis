# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:31 2017

@author: mondon
"""

#! /usr/bin/env python
import sncosmo
import builtins_SNF
import os
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import read_data
import numpy as np
import copy
import read_data as RD
#from decimal import Decimal

t_min = -15
t_max = 45

t_min_sug = -12
t_max_sug = 42

wl_min_sal = 3000.
wl_max_sal = 7000.

def fit_salt2(meta):

#    outfile = open('../sugar_analysis_data/results/res_salt2_SNF.txt', 'w')
#    outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor cov_m_s cov_m_c cov_s_c  tmax dtmax x0 dx0')

    list_SN = ['SNF20080323-009']
   
    
    
    fitfail=[]   

#    for sn_name in meta.keys():
    for sn_name in list_SN:
#        if 'lc-SDSS' == filename[:7] or 'lc-sn'==filename[:5]:
#        if 'lc-SDSS' == filename[:7]:   
        
        
        print sn_name
        data, zhl, zcmb, mwebv = RD.read_meta_SNF(meta,sn_name,filters=['BSNf','VSNf','RSNf'])
        source = sncosmo.get_source('salt2', version='2.4')
        source.EBV_snfit = mwebv    
        source.z_snfit = zhl
        source.Rv_snfit = 3.1

        dust = sncosmo.CCM89Dust()
        model = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])
#        model = sncosmo.Model(source=source)

        model.set(mwebv=mwebv)
        model.set(z=zhl)
              
#        try:     
        #initial iteration x1 fix
        model.set(x1=0.) #with 0.01 the  fit work for all SDSS and nearby

        try:    
            print 'initialisation'
            res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'c'], modelcov = False, apply_ZPERR=False)
            print 'first iteration'
            chi2 = res.chisq
            print chi2
            chi2p = chi2*2
            m=0

            while chi2 < chi2p:


                if m > 0:
                    resp = res
                    fitted_modelp = fitted_model
                m += 1
                t_peak = fitted_model.parameters[1]
                #print t_peak,fitted_model.parameters[4]
    
                t1 = t_peak + t_min*(1 + model.get('z'))
                t2 = t_peak + t_max*(1 + model.get('z'))
                            
                A=[]
                data_new = copy.deepcopy(data)
                for i in range(len(data_new)):                    
                    if data_new[i][0] <= t1 or data_new[i][0] >= t2:
                        A.append(i)
                A=np.array(A)
                for i in range(len(A)):
                    #print  'We excluded the point %7.3f because it does not belong to the time interval [%7.2f,%7.2f]' % (data_new[A[i]][0],t1,t2)
                    data_new.remove_row(A[i])
                    A-=1   
               
                res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0', 'x0', 'x1', 'c'], modelcov = True, apply_ZPERR=False) 
                chi2p = chi2
    
                chi2 = res.chisq
                print chi2p, chi2
               

            
                
            #final results
            res = resp
            fitted_model = fitted_modelp
            print (res.chisq)
            
            #Calculation of mb
            mb = mB_determination(res)
            print mb

    #                
    
    #        
            #Calculation of mb uncertainty and covariance between mb and other parameters
            dmbfit, cov_mb_c, cov_mb_x1, cov_x1_c, cov_mb_x0 = mb_uncertainty(res)
    #        print dmbfit, cov_mb_c, cov_mb_x1, cov_x1_c, cov_mb_x0 
    #        
#        
#
#        print res, '\nNumber of iterations: ', m
            sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
            plt.show()

#            outfile.write('\n')        
#            outfile.write('%s 999 %f 999 %f %f %f %f %f %f %e %e %e %f %f %f %f' %(sn_name, res.parameters[0], mb, dmbfit, res.parameters[3], res.errors['x1'], res.parameters[4], res.errors['c'], cov_mb_x1, cov_mb_c, cov_x1_c, res.parameters[2], res.errors['x0'], res.parameters[1], res.errors['t0']))          

#
        except:
            fitfail.append(sn_name)
            print 'Error: fit fail for: ',sn_name          
    
            
    print fitfail       
#    outfile.close()
    return res

#    
    
def fit_sugar():

    #outfile = open('res_sugar.txt', 'w')
    for filename in os.listdir('/home/maria/Dropbox/Science/Supernovae/JLA_SALT2/JLA_fit/jla_data'):
    #list_SDSS = ['lc-SDSS18468.list']
    #list_SDSS = ['b_lc-SDSS14318.list','b_lc-SDSS17274.list','b_lc-SDSS20470.list','b_lc-SDSS21042.list','b_lc-SDSS16402.list','b_lc-SDSS11300.list','b_lc-SDSS16206.list','b_lc-SDSS20575.list','b_lc-SDSS16281.list','b_lc-SDSS18617.list','b_lc-SDSS17220.list']    
    #for filename in list_SDSS:
        if 'lc-SDSS' == filename[:7]:
            sn_name = filename
            print sn_name
            source = sncosmo.get_source('sugar')    
            dust = sncosmo.CCM89Dust()
            model = sncosmo.Model(source=source, effects=[dust], effect_names=['mw'], effect_frames=['obs'])
            
            head = read_lc_jla(sn_name)
            model.set(mwebv=head['@MWEBV'])
            model.set(z=head['@Z_HELIO'])
    
            data = sncosmo.read_lc('/home/maria/Dropbox/Science/Supernovae/sncosmo/test_jla/' + sn_name)
    
            print head['@DayMax']
            res, fitted_model = sncosmo.fit_lc(data, model, ['t0','q1','q2', 'q3', 'A', 'Mgr'], bounds={'t0':(head['@DayMax']-3,head['@DayMax']+3)},guess_t0=True,phase_range=(-12,42))

            print res
            sncosmo.plot_lc(data, model=fitted_model, errors=res.errors)
            plt.show()

            #outfile.write('%s \n' %(sn_name)) 
            #outfile.write('%s \n' %(res)) 
        
    #outfile.close()    


    


def mB_determination(res):
    
    scale_factor = 10**-12
    
    #interpolation of TB and Trest
    filt2 = np.genfromtxt('../sncosmo_jla/jla_data/Instruments/SNLS3-Landolt-model/sb-shifted.dat')
    wlen = filt2[:,0]
    tran = filt2[:,1]
    splB = Spline1d(wlen, tran, k=1,ext = 1)



    #interpolation of ref spectrum
    data = np.genfromtxt('../sncosmo_jla/jla_data/MagSys/bd_17d4708_stisnic_002.ascii')
    dispersion = data[:,0]
    flux_density = data[:,1]
    splref = Spline1d(dispersion, flux_density, k=1,ext = 1)

  
    #interpolation of the spectrum model
    template_0 = np.genfromtxt('../sncosmo_jla/salt2-4/salt2_template_0.dat')    
    template_1 = np.genfromtxt('../sncosmo_jla/salt2-4/salt2_template_1.dat')
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

    return mb

def color_law_salt2(wl):
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

           
def mb_uncertainty(res):
    
    h = 10**-9
    
    #build mb derivative for all parameters    
    resx0 = copy.deepcopy(res)  
    resx0.parameters[2] = resx0.parameters[2] + h
    dmb_dx0 = (mB_determination(resx0)- mB_determination(res))/h

    resx1 = copy.deepcopy(res)  
    resx1.parameters[3] = resx0.parameters[3] + h
    dmb_dx1 = (mB_determination(resx1)- mB_determination(res))/h

    resc = copy.deepcopy(res)  
    resc.parameters[4] = resx0.parameters[4] + h
    dmb_dc = (mB_determination(resc)- mB_determination(res))/h
    
    #build vectors for all mb derivative
    vect = np.array([dmb_dx0, dmb_dx1, dmb_dc])
    
    #build the covariance matrix for salt2 parameters

    mat = np.delete(res.covariance, (0), axis=0)
    mat = mat.T
    mat = np.delete(mat, (0), axis=0)
    mat = mat.T
#    mat = res.covariance
#    print mat
#    print len(mat)
#    print mat
    
    dmb = np.sqrt(np.dot(np.dot(vect.T, mat), vect))
    cov_mb_c = mat[2,2] * dmb_dc + mat[2,1] * dmb_dx1 + mat[2,0] * dmb_dx0
    cov_mb_x1 = mat[1,1] * dmb_dx1 + mat[1,0] * dmb_dx0 + mat[1,2] * dmb_dc
    cov_mb_x0 = mat[0,0] * dmb_dx0 +  mat[0,1] * dmb_dx1 + mat[0,2] * dmb_dc
    cov_x1_c = mat[2,1]
    return dmb, cov_mb_c, cov_mb_x1, cov_x1_c, cov_mb_x0