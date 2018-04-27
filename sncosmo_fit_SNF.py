# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 13:24:31 2017

@author: mondon
"""

#! /usr/bin/env python
import sncosmo
import builtins_SNF as Build_SNF
import os
from matplotlib import pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.gridspec as gridspec
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import numpy as np
import copy
import cPickle
import read_data as RD
#from decimal import Decimal

t_min = -15
t_max = 45

t_min_sug = -12
t_max_sug = 48

wl_min_sal = 3000.
wl_max_sal = 7000.

def Light_curve_fit(filters=['BSNf','VSNf','RSNf'],errorscale=True, model_used='salt2', modelcov=True, width=10, write_results=False):
    
    if write_results :
        if model_used=='salt2':
            if modelcov:
                outfile = open('../sugar_analysis_data/results/res_salt2_SNF_'+str(width)+'.txt', 'w')
            else:
                outfile = open('../sugar_analysis_data/results/res_salt2_SNF_'+str(width)+'_nomodelcov.txt', 'w')
            outfile.write('#name zcmb zhel dz mb dmb x1 dx1 color dcolor cov_m_s cov_m_c cov_s_c x0 dx0 tmax dtmax chi2')
        elif model_used == 'sugar':
            outfile = open('../sugar_analysis_data/results/res_sugar_SNF.txt', 'w')
            outfile.write('#name zcmb zhel dz mgrey dmgrey q1 dq1 q2 dq2 q3 dq3 A dA cov_mgrey_q1 cov_mgrey_q2 cov_mgrey_q3 cov_mgrey_A cov_q1_q2 cov_q1_q3 cov_q1_A cov_q2_q3 cov_q2_A cov_q3_A tmax dtmax chi2')
        else:
             raise ValueError('ERROR: model name has to be salt2 or sugar')
            
#    list_SN = ['SNF20070531-011']
    list_SN = ['SN2008ec','SNF20080522-011','PTF11qzq']

#    list_SN = ['SN2008ec']
#    list_SN = ['SNF20080510-001','SN2005ir','SNF20080919-000','SNF20080522-011', 'SNF20070417-002','SNF20070424-003','SNF20080610-000','SN2007nq','SN2008ec']
    meta = cPickle.load(open('../sugar_analysis_data/META-CABALLO2.pkl'))
    
    
    fitfail=[]   

#    for sn_name in meta.keys():

    for sn_name in list_SN :
        if meta[sn_name]['idr.subset'] != 'bad' and meta[sn_name]['idr.subset'] != 'auxiliary' and sn_name != 'SN2011by' and sn_name != 'SNNGC4424' and sn_name != 'SN2012fr' :
  
        
        
            print sn_name
            data, zhl, mwebv, daymax = RD.read_meta_SNF(meta, sn_name, filters=filters, errorscale=errorscale, width=width)

            if model_used == 'salt2':          
                source = sncosmo.get_source('salt2', version='2.4')
            elif model_used == 'sugar':
                source = sncosmo.get_source('sugar')
            else:
                raise ValueError('ERROR: model name has to be salt2 or sugar')
#            source.EBV_snfit = mwebv    
#            source.z_snfit = zhl
#            source.Rv_snfit = 3.1
    
            dust = sncosmo.CCM89Dust()
            model = sncosmo.Model(source=source, effects=[dust],effect_names=['mw'], effect_frames=['obs'])    
            model.set(mwebv=mwebv)
            model.set(z=zhl)
                  
            try:     
                    #initial iteration x1 fix
        
            
                print 'initialisation'
                
                
                if model_used =='salt2':      
                    if sn_name == 'SNF20080717-000':
                        model.set(x1=3.)
                    else:
                        model.set(x1=0.3)
                    
                    res, fitted_model = sncosmo.fit_lc(data, model, ['t0','x0', 'c'], modelcov = False)
                elif model_used == 'sugar':
    

                    model.set(q1=0.0)
                    model.set(q2=0.0)
                    model.set(q3=0.0)     
                    model.set(t0=daymax)
    #                    print daymax
                    res, fitted_model = sncosmo.fit_lc(data, model, ['A','Mgr'], modelcov = False)
    #                    print res.chisq
    #                    model.set(q1=res.parameters[2])
    #                    model.set(q2=res.parameters[3])
    #                    model.set(q3=res.parameters[4])
    #                model.set(A=res.parameters[5])
    #                model.set(Mgr=res.parameters[6])
    #                res, fitted_model = sncosmo.fit_lc(data, model, ['q1', 'q2', 'q3', 'A', 'Mgr'], modelcov  = False)
    #                    res, fitted_model = sncosmo.fit_lc(data, model, ['t0', 'A', 'Mgr'], modelcov  = False,phase_range=(-12,42))
                else:
                    raise ValueError('ERROR: model name has to be salt2 or sugar')
                    
                print 'first iteration'
                chi2 = res.chisq
                print chi2
                chi2p = chi2*2
                m=0
            
                while chi2 < chi2p and m < 10:
        
        
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
        
                    if model_used =='salt2':          
                        res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0', 'x0', 'x1', 'c'], modelcov = modelcov)
                    elif model_used == 'sugar':
                        print "WARNING: Sugar don't have model covariance for the moment"
                        modelcov = False
                        model.set(q1=res.parameters[2])
                        model.set(q2=res.parameters[3])
                        model.set(q3=res.parameters[4])
                        model.set(A=res.parameters[5])
                        model.set(Mgr=res.parameters[6])
    #                        print res.parameters[6]
                        model.set(t0=res.parameters[1])                        
                        res, fitted_model = sncosmo.fit_lc(data_new, model, ['t0','q1', 'q2', 'q3', 'A', 'Mgr'], modelcov  = modelcov)
                        print res.parameters[1]
    
#                        sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
                        
    #                        print res.chisq
        #                    print res
    
                    else:
                        raise ValueError('ERROR: model name has to be salt2 or sugar')                   
    #                    print m, chi2p, chi2
                    chi2p = chi2
                
                    chi2 = res.chisq
                    if m == 1 and chi2p <= chi2:
                        raise ValueError('ERROR: at the first iteration chi2p have to be inferior to chi2')
    
                    print chi2p, chi2
                    
        
                
                    
                #final results
                res = resp
                fitted_model = fitted_modelp
                print (res.chisq)
    #                print res
                
                #Calculation of mb
                
                if model_used =='salt2':          
                    mb = mB_determination(res)
                    print mb
                    #Calculation of mb uncertainty and covariance between mb and other parameters
                    dmbfit, cov_mb_c, cov_mb_x1, cov_x1_c, cov_mb_x0 = mb_uncertainty(res)
                elif model_used == 'sugar':
                    print 'For the moment no mb for sugar use Mgrey'
#                    sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
#                    plt.savefig('../sugar_analysis_data/results/'+sn_name+'_sugar_fit.pdf')
#                    plt.show()
        
                    cov_mgrey_q1, cov_mgrey_q2, cov_mgrey_q3, cov_mgrey_A, cov_q1_q2, cov_q1_q3, cov_q1_A, cov_q2_q3, cov_q2_A, cov_q3_A = sugar_covariance(res)
                else:
                    raise ValueError('ERROR: model name has to be salt2 or sugar')      
                    
                if write_results :
                    outfile.write('\n')        
                    if model_used =='salt2':          
                        outfile.write('%s 999 %f 999 %f %f %f %f %f %f %e %e %e %f %f %f %f %f' %(sn_name, res.parameters[0], mb, dmbfit, res.parameters[3], res.errors['x1'], res.parameters[4], res.errors['c'], cov_mb_x1, cov_mb_c, cov_x1_c, res.parameters[2], res.errors['x0'], res.parameters[1], res.errors['t0'], chi2))          
                    elif model_used == 'sugar':  
#                        outfile.write('%s 999 %f 999 %f %f %f %f %f %f %f %f %f %f %e %e %e %e %e %e %e %e %e %e %f' %(sn_name, res.parameters[0], res.parameters[6], res.errors['Mgr'], res.parameters[2], res.errors['q1'], res.parameters[3], res.errors['q2'], res.parameters[4], res.errors['q3'],res.parameters[5], res.errors['A'], cov_mgrey_q1, cov_mgrey_q2, cov_mgrey_q3, cov_mgrey_A, cov_q1_q2, cov_q1_q3, cov_q1_A, cov_q2_q3, cov_q2_A, cov_q3_A, chi2))
                        outfile.write('%s 999 %f 999 %f %f %f %f %f %f %f %f %f %f %e %e %e %e %e %e %e %e %e %e %f %f %f' %(sn_name, res.parameters[0], res.parameters[6], res.errors['Mgr'], res.parameters[2], res.errors['q1'], res.parameters[3], res.errors['q2'], res.parameters[4], res.errors['q3'],res.parameters[5], res.errors['A'], cov_mgrey_q1, cov_mgrey_q2, cov_mgrey_q3, cov_mgrey_A, cov_q1_q2, cov_q1_q3, cov_q1_A, cov_q2_q3, cov_q2_A, cov_q3_A, res.parameters[1], res.errors['t0'], chi2))
    
    #                        outfile.write('%s 999 %f 999 %f %f %f %f %f %f %f %f %f %f 999 999 999 999 999 999 999 999 999 999 %f %f %f' %(sn_name, res.parameters[0], res.parameters[6], res.errors['Mgr'], res.parameters[2], res.errors['q1'], res.parameters[3], res.errors['q2'], res.parameters[4], res.errors['q3'],res.parameters[5], res.errors['A'], res.parameters[1], 0, chi2))
                    else:
                        raise ValueError('ERROR: model name has to be salt2 or sugar')      
#
            except:
                fitfail.append(sn_name)
                print 'Error: fit fail for: ',sn_name          
    
            
    print fitfail
    print len(fitfail)
    SUGAR_parameter_pkl = '../sugar/sugar/data_output/data_output/sugar_parameters.pkl'
    dico = cPickle.load(open(SUGAR_parameter_pkl)) 
    for sn_fail in fitfail:
        if sn_fail in dico.keys():
            print sn_fail
        
    if write_results :           
        outfile.close()
    return res , data_new, model , data

def mB_determination(res):
    """
    determine mb restframe 
    """
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
#    fits = pyfits.open('../sugar_analysis_data/Vega.fits')
#    fit = fits[1]
#    wl = np.zeros(len(fit.data))
#    flux = np.zeros(len(fit.data))
#    for i in range(len(fit.data)):
#         wl[i] = fit.data[i][0]
#         flux[i] = fit.data[i][1]

    splref = Spline1d(dispersion, flux_density, k=1,ext = 1)

  
    #interpolation of the spectrum flux
    template_0 = np.genfromtxt('../sncosmo_jla/jla_data/salt2-4/salt2_template_0.dat')    
    template_1 = np.genfromtxt('../sncosmo_jla/jla_data/salt2-4/salt2_template_1.dat')
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

def sugar_covariance(res):
    '''
    '''
    
#    cov_mgrey_q1 = res.covariance[5,1]
#    cov_mgrey_q2 = res.covariance[5,2]
#    cov_mgrey_q3 = res.covariance[5,3]
#    cov_mgrey_A = res.covariance[5,4]
#    cov_q1_q2 = res.covariance[1,2]
#    cov_q1_q3 = res.covariance[1,3]
#    cov_q1_A = res.covariance[1,4]
#    cov_q2_q3 = res.covariance[2,3]
#    cov_q2_A = res.covariance[2,4]
#    cov_q3_A = res.covariance[3,4]  
    
    cov_mgrey_q1 = res.covariance[4,0]
    cov_mgrey_q2 = res.covariance[4,1]
    cov_mgrey_q3 = res.covariance[4,2]
    cov_mgrey_A = res.covariance[4,3]
    cov_q1_q2 = res.covariance[0,1]
    cov_q1_q3 = res.covariance[0,2]
    cov_q1_A = res.covariance[0,3]
    cov_q2_q3 = res.covariance[1,2]
    cov_q2_A = res.covariance[1,3]
    cov_q3_A = res.covariance[2,3]      
    
    return cov_mgrey_q1, cov_mgrey_q2, cov_mgrey_q3, cov_mgrey_A, cov_q1_q2, cov_q1_q3, cov_q1_A, cov_q2_q3, cov_q2_A, cov_q3_A

def mb_uncertainty(res):
    """
    determine mb uncertainty
    """
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

    dmb = np.sqrt(np.dot(np.dot(vect.T, mat), vect))
    cov_mb_c = mat[2,2] * dmb_dc + mat[2,1] * dmb_dx1 + mat[2,0] * dmb_dx0
    cov_mb_x1 = mat[1,1] * dmb_dx1 + mat[1,0] * dmb_dx0 + mat[1,2] * dmb_dc
    cov_mb_x0 = mat[0,0] * dmb_dx0 +  mat[0,1] * dmb_dx1 + mat[0,2] * dmb_dc
    cov_x1_c = mat[2,1]
    return dmb, cov_mb_c, cov_mb_x1, cov_x1_c, cov_mb_x0

if __name__=="__main__":
    
    width = 10
    model_used = 'sugar'
#    model_used = 'salt2'
    write_results = False
    modelcov = True
    try:
        Build_SNF.register_SNf_bands_width(width=width)
        Build_SNF.mag_sys_SNF_width(width=width)
        Build_SNF.register_SUGAR()    
    except:
        print 'Filters and mag sys already registred' 
    
    res, data_new, model,data = Light_curve_fit(model_used=model_used,width=width, write_results=write_results, modelcov=modelcov)