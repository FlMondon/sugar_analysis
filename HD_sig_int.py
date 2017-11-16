# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 12:20:25 2017

@author: mondon
"""

import cPickle
import numpy as np
import pylab as P
import sugar
from scipy import integrate
from numpy.linalg import inv
import iminuit as minuit
from scipy import optimize,integrate




###########
#constants
###########

clight = 299792.458
H0 = 0.000070
###############



def read_sugar_salt2_par():
    
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
    
    
    zcmb = zcmb[Filtre]
    zhl = zhl[Filtre]
    zerr = zerr[Filtre]
    cov_x = cov_x[Filtre]
    grey = grey[Filtre]
    q1 = q1[Filtre]
    q2 = q2[Filtre]
    q3 = q3[Filtre]
    av = av[Filtre]
    
      
    X = np.array([grey,q1,q2,q3,av]).T
    Y = np.array([mb,x1,c]).T
    print len(Y)
    return X, Y, cov_x, cov_y,  zhl, zcmb, zerr

def comp_rms(residuals, dof, err=True, variance=None):
    """
    Compute the RMS or WRMS of a given distribution.

    :param 1D-array residuals: the residuals of the fit.
    :param int dof: the number of degree of freedom of the fit.
    :param bool err: return the error on the RMS (WRMS) if set to True.
    :param 1D-aray variance: variance of each point. If given,
                             return the weighted RMS (WRMS).

    :return: rms or rms, rms_err
    """
    if variance is None:                # RMS
        rms = float(np.sqrt(np.sum(residuals**2)/dof))
        rms_err = float(rms / np.sqrt(2*dof))
    else:                               # Weighted RMS
        assert len(residuals) == len(variance)
        rms = float(np.sqrt(np.sum((residuals**2)/variance) / np.sum(1./variance)))
        #rms_err = float(N.sqrt(1./N.sum(1./variance)))
        rms_err = np.sqrt(2.*len(residuals)) / (2*np.sum(1./variance)*rms)

    if err:
        return rms, rms_err
    else:
        return rms
        
class Hubble_fit():   
    
    def int_cosmo(self, z, Omega_M=0.3):     
        return 1./np.sqrt(Omega_M*(1+z)**3+(1.-Omega_M))
        
    def luminosity_distance(self):
        
        if type(self.zcmb)==np.ndarray:
            integr = np.zeros_like(self.zcmb)
            for i in range(len(self.zcmb)):
                integr[i]=integrate.quad(self.int_cosmo, 0, self.zcmb[i])[0]
        else:
            integr = integrate.quad(self.int_cosmo, 0, self.zcmb)[0]
    
        return (1+self.zhl)*(clight/H0)*integr
 
    def distance_modulus_th(self):      
        return 5.*np.log(self.luminosity_distance())/np.log(10.)-5.
        


class Hubble_fit_salt2(Hubble_fit):
    
    def __init__(self, Y, cov_y, zhl, zcmb, zerr):
        self.salt2_par = Y
        self .cov = cov_y
        self.zcmb = zcmb
        self.zhl = zhl
        self.zerr = zerr
        self.dmz = 5/np.log(10) * np.sqrt(self.zerr**2 + 0.001**2) / self.zcmb
        self.dof = len(Y)-3        

    def distance_modulus(self, alpha, beta, Mb):
        return self.salt2_par[:,0] - Mb + alpha*self.salt2_par[:,1] - beta*self.salt2_par[:,2] 

    def chi2(self, alpha, beta, Mb, sig_int):


        Cmu = np.zeros_like(self.cov[::3,::3])
        for i, coef1 in enumerate([1., alpha, -beta]):
            for j, coef2 in enumerate([1., alpha, -beta]):
                Cmu += (coef1 * coef2) * self.cov[i::3,j::3]               
        Cmu[np.diag_indices_from(Cmu)] += sig_int**2 + self.dmz**2 
        C = inv(Cmu)
        L = self.distance_modulus(alpha, beta, Mb) - self.distance_modulus_th()
        self.residuals = L
        self.var = np.diag(Cmu)
        return P.dot(L,P.dot(C,L))
        
    def _compute_dispertion(self, sig_int):
        sig = optimize.fmin(self._disp_function, sig_int)[0]
        return sig
             
    def _disp_function(self,d):
        return abs((self.chi2(self.Params['alpha'], self.Params['beta'], self.Params['Mb'], d)/self.dof)-1.)
               
         
    def fit_sigma_int(self):
        
        sig_int = 0.001
        
        Find_param = minuit.Minuit(self.chi2, alpha=0.15,beta=2.9, Mb=-19.1, sig_int=sig_int, fix_sig_int= True)
        Find_param.migrad()
        self.Params = Find_param.values
        self.Params_Covariance = Find_param.covariance
        

        calls=0
        if abs((self.chi2(self.Params['alpha'], self.Params['beta'], self.Params['Mb'], sig_int)/(self.dof))-1.)>0.1:
            while abs((self.chi2(self.Params['alpha'], self.Params['beta'], self.Params['Mb'], sig_int)/(self.dof))-1.) > 0.001:
   
                if calls<100:
                    print 'je cherche de la dispersion pour la %i eme fois'%(calls+1)
                    sig_int = self._compute_dispertion(sig_int)                  
                    Find_param = minuit.Minuit(self.chi2, alpha=self.Params['alpha'], beta=self.Params['beta'], Mb=self.Params['Mb'], sig_int = sig_int, fix_sig_int= True)
                
                    Find_param.migrad()
                    self.Params = Find_param.values
                    print self.Params
                    calls+=1

                else:
                    print 'error : calls limit are exceeded'
                    break
        
        self.wrms, self.wrms_err = comp_rms(self.residuals, self.dof, err=True, variance=self.var)
        return Find_param, sig_int        
        
class Hubble_fit_sugar(Hubble_fit):

    def __init__(self, X, cov_x, zhl, zcmb, zerr, qi=True):
        self.sugar_par = X
        self.cov = cov_x
        self.zcmb = zcmb
        self.zhl = zhl
        self.zerr = zerr
        self.dmz = 5/np.log(10) * np.sqrt(self.zerr**2 + 0.001**2) / self.zcmb       
        self.qi = qi
        
        if self.qi == False:
            self.dof = len(X) - 1
        else:
            self.dof = len(X)- 5
    
    def distance_modulus_grey(self, cst):      
        return self.sugar_par[:,0] - cst

    def distance_modulus_corr(self, cst, alpha1, alpha2, alpha3, beta):      
        return self.sugar_par[:,0] - cst - alpha1*self.sugar_par[:,1] - alpha2*self.sugar_par[:,2] - alpha3*self.sugar_par[:,3] - beta*self.sugar_par[:,4]
        

    
    def build_cov_mat(self):
        
        if self.qi == False:
            cov_mat = np.zeros([len(self.sugar_par), len(self.sugar_par)])
            cv = self.cov
    
            for i in range (len(self.sugar_par)):
                cov_mat[i,i] = cv[i,0,0]
            self.cov_mat = cov_mat
        else: 
            cov_mat = np.zeros([len(self.sugar_par)*5, len(self.sugar_par)*5])
            cv = self.cov
            for i in range (len(self.sugar_par)):
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
                self.cov_mat  = cov_mat 
                
        
    
    

    def chi2(self, cst, alpha1, alpha2, alpha3, beta,sig_int):
        

        if self.qi == False:
            self.build_cov_mat()
            Cmu =  self.cov_mat
            Cmu[np.diag_indices_from(Cmu)] += sig_int**2 + self.dmz**2 
            
            L = self.distance_modulus_grey(cst) #- self.distance_modulus_th()

        else:
            self.build_cov_mat()
            Cmu = np.zeros_like(self.cov_mat[::5,::5])
            for i, coef1 in enumerate([1., -alpha1,-alpha2,-alpha3, -beta]):
                for j, coef2 in enumerate([1., -alpha1,-alpha2,-alpha3, -beta]):
                    Cmu += (coef1 * coef2) * self.cov_mat[i::5,j::5]
            Cmu[np.diag_indices_from(Cmu)] += sig_int**2 + self.dmz**2 
            L = self.distance_modulus_corr(cst, alpha1, alpha2, alpha3, beta) #- self.distance_modulus_th()
            
        self.Cmu = Cmu
        C = inv(Cmu)
        self.residuals = L
        self.var = np.diag(Cmu)            
        return P.dot(L.T,P.dot(C,L))        
        








    def _compute_dispertion(self, sig_int):
        sig = optimize.fmin(self._disp_function, sig_int)[0]
        print sig
        return sig
             
    def _disp_function(self,d):
        print d
        return abs((self.chi2(self.Params['cst'], self.Params['alpha1'], self.Params['alpha2'], self.Params['alpha3'], self.Params['beta'] ,d)/self.dof)-1.)        
        
        
    def fit_sigma_int(self):
        sig_int = 0.001
        
        if self.qi == False: 
            Find_param = minuit.Minuit(self.chi2, cst=0.01,alpha1=0., alpha2=0., alpha3=0., beta=0.,sig_int=sig_int, fix_sig_int= True, fix_alpha1=True, fix_alpha2=True, fix_alpha3=True, fix_beta=True)
            print sig_int
            Find_param.migrad()
            self.Params = Find_param.values
            self.Params_Covariance = Find_param.covariance
            
    
            calls=0
            if abs((self.chi2(self.Params['cst'], 0., 0., 0., 0.,sig_int)/(self.dof))-1.)>0.1:
                while abs((self.chi2(self.Params['cst'],0., 0., 0., 0.,sig_int)/(self.dof))-1.) > 0.001:
       
                    if calls<100:
                        print 'je cherche de la dispersion pour la %i eme fois'%(calls+1)
                        sig_int = self._compute_dispertion(sig_int)
                        print sig_int
                        Find_param = minuit.Minuit(self.chi2, cst=self.Params['cst'], alpha1=0., alpha2=0., alpha3=0., beta=0., sig_int = sig_int, fix_sig_int=True, fix_alpha1=True, fix_alpha2=True, fix_alpha3=True, fix_beta=True)             
                        Find_param.migrad()
                        self.Params = Find_param.values
                        print self.Params
                        calls+=1
                        
    
                    else:
                        print 'error : calls limit are exceeded'
                        break
        else:
            Find_param = minuit.Minuit(self.chi2, cst=0.01 ,alpha1=0.001, alpha2=0.001, alpha3=0.001, beta=0.001, sig_int=sig_int, fix_sig_int= True)
            print sig_int
            Find_param.migrad()
            self.Params = Find_param.values
            self.Params_Covariance = Find_param.covariance
            
    
            calls=0
            if abs((self.chi2(self.Params['cst'], self.Params['alpha1'], self.Params['alpha2'],self.Params['alpha3'], self.Params['beta'],sig_int)/(self.dof))-1.)>0.1:
                while abs((self.chi2(self.Params['cst'], self.Params['alpha1'], self.Params['alpha2'],self.Params['alpha3'], self.Params['beta'], sig_int)/(self.dof))-1.) > 0.001:
       
                    if calls<100:
                        print 'je cherche de la dispersion pour la %i eme fois'%(calls+1)
                        sig_int = self._compute_dispertion(sig_int)
                        print sig_int
                        Find_param = minuit.Minuit(self.chi2, cst=self.Params['cst'], alpha1=self.Params['alpha1'], alpha2=self.Params['alpha2'], alpha3=self.Params['alpha3'], beta=self.Params['beta'], sig_int = sig_int, fix_sig_int= True)             
                        Find_param.migrad()
                        self.Params = Find_param.values
                        print self.Params
                        calls+=1
                        
    
                    else:
                        print 'error : calls limit are exceeded'
                        break
            
        self.wrms, self.wrms_err = comp_rms(self.residuals, self.dof, err=True, variance=self.var)        
        return Find_param, sig_int
        