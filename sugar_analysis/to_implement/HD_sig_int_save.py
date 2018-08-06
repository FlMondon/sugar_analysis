# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 12:20:25 2017

@author: mondon
"""

import cPickle
import numpy as np
import pylab as P
import sugar
from matplotlib import rc, rcParams
from scipy import integrate
from numpy.linalg import inv
import iminuit as minuit
from scipy import optimize,integrate
from matplotlib import pyplot as plt
import copy



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
        
class Hubble_fit():   
    
    def __init__(self, X, Y, cov_x, cov_y, zhl, zcmb, zerr, qi=True):
        self.sugar_par = X
        self.salt2_par = Y
        self.cov_sugar = cov_x
        self .cov_salt2 = cov_y
        self.zcmb = zcmb
        self.zhl = zhl
        self.zerr = zerr
        self.dmz = 5/np.log(10) * np.sqrt(self.zerr**2 + 0.001**2) / self.zcmb
        self.dof_salt2 = len(Y)-3        
        self.qi = qi
        print self.cov_sugar
        if self.qi == False:
            self.dof_sugar = len(X) - 1
        else:
            self.dof_sugar = len(X)- 5
            
            
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
 
    def distance_modulus(self):      
        return 5.*np.log(self.luminosity_distance())/np.log(10.)-5.
        
    def distance_modulus_sugar_grey(self, cst):      
        return self.sugar_par[:,0] - cst

    def distance_modulus_sugar_corr(self, cst, alpha1, alpha2, alpha3, beta):      
        return self.sugar_par[:,0] - cst - alpha1*self.sugar_par[:,1] - alpha2*self.sugar_par[:,2] - alpha3*self.sugar_par[:,3] - beta*self.sugar_par[:,4]
    
    def distance_modulus_sugar_corr_uncertainty(self, alpha1, alpha2, alpha3, beta, dMgr, dq1, dq2, dq3, dA):
        return np.sqrt(dMgr**2 + (alpha1*dq1)**2 + (alpha2*dq2)**2 + (alpha3*dq3)**2 + (beta*dA)**2)
    
    def distance_modulus_salt2(self, alpha, beta, Mb):
        return self.salt2_par[:,0] - Mb + alpha*self.salt2_par[:,1] - beta*self.salt2_par[:,2] 
    
    def build_cov_mat_sugar(self):
        
        if self.qi == False:
            cov_mat_sugar = np.zeros([len(self.sugar_par), len(self.sugar_par)])
            cv = self.cov_sugar
    
            for i in range (len(self.sugar_par)):
                cov_mat_sugar[i,i] = cv[i,0,0]
            self.cov_mat_sugar = cov_mat_sugar
        else: 
            cov_mat_sugar = np.zeros([len(self.sugar_par)*5, len(self.sugar_par)*5])
            cv = self.cov_sugar
            for i in range (len(self.sugar_par)):
                cov_mat_sugar[i*5,i*5] = cv[i,0,0]
                cov_mat_sugar[i*5 +1,i*5] = cv[i,1,0]
                cov_mat_sugar[i*5 +2,i*5] = cv[i,2,0]
                cov_mat_sugar[i*5 +3,i*5] = cv[i,3,0]
                cov_mat_sugar[i*5 +4,i*5] = cv[i,4,0]
                cov_mat_sugar[i*5,i*5 +1] = cv[i,0,1]
                cov_mat_sugar[i*5 +1,i*5 +1] = cv[i,1,1]
                cov_mat_sugar[i*5 +2,i*5 +1] = cv[i,2,1]
                cov_mat_sugar[i*5 +3,i*5 +1] = cv[i,3,1]
                cov_mat_sugar[i*5 +4,i*5 +1] = cv[i,4,1]
                cov_mat_sugar[i*5 ,i*5 +2] = cv[i,0,2]
                cov_mat_sugar[i*5 +1,i*5 +2] = cv[i,1,2]
                cov_mat_sugar[i*5 +2,i*5 +2] = cv[i,2,2]
                cov_mat_sugar[i*5 +3,i*5 +2] = cv[i,3,2]
                cov_mat_sugar[i*5 +4,i*5 +2] = cv[i,4,2]
                cov_mat_sugar[i*5,i*5 +3] = cv[i,0,3]
                cov_mat_sugar[i*5 +1,i*5 +3] = cv[i,1,3]
                cov_mat_sugar[i*5 +2,i*5 +3] = cv[i,2,3]
                cov_mat_sugar[i*5 +3,i*5 +3] = cv[i,3,3]
                cov_mat_sugar[i*5 +4,i*5 +3] = cv[i,4,3]
                cov_mat_sugar[i*5,i*5 +4] = cv[i,0,4]
                cov_mat_sugar[i*5 +1,i*5 +4] = cv[i,1,4]
                cov_mat_sugar[i*5 +2,i*5 +4] = cv[i,2,4]
                cov_mat_sugar[i*5 +3,i*5 +4] = cv[i,3,4]
                cov_mat_sugar[i*5 +4,i*5 +4] = cv[i,4,4]
                self.cov_mat_sugar = cov_mat_sugar
                
        
    
    

    def chi2_sugar(self, cst, alpha1, alpha2, alpha3, beta,sig_int):
        

        if self.qi == False:
#            self.build_cov_mat_sugar()
#            self.cov_mat_sugar = self.cov_sugar
            Cmu =  copy.deepcopy(self.cov_sugar)
            Cmu[np.diag_indices_from(Cmu)] += sig_int**2 + self.dmz**2 
            
            L = self.distance_modulus_sugar_grey(cst) - self.distance_modulus()

        else:
#            self.build_cov_mat_sugar()
#            self.cov_mat_sugar = self.cov_sugar
            Cmu = np.zeros_like(self.cov_sugar[::5,::5])
            for i, coef1 in enumerate([1., -alpha1,-alpha2,-alpha3, -beta]):
                for j, coef2 in enumerate([1., -alpha1,-alpha2,-alpha3, -beta]):
                    Cmu += (coef1 * coef2) * self.cov_sugar[i::5,j::5]
            Cmu[np.diag_indices_from(Cmu)] += sig_int**2 + self.dmz**2 
            L = self.distance_modulus_sugar_corr(cst, alpha1, alpha2, alpha3, beta) - self.distance_modulus()
            
        self.Cmu = Cmu
        C = inv(Cmu)
        self.residuals = L
        self.var = np.diag(Cmu)            
        return P.dot(L.T,P.dot(C,L))        
        
    def chi2_salt2(self, alpha, beta, Mb, sig_int):


        Cmu = np.zeros_like(self.cov_salt2[::3,::3])
        for i, coef1 in enumerate([1., alpha, -beta]):
            for j, coef2 in enumerate([1., alpha, -beta]):
                Cmu += (coef1 * coef2) * self.cov_salt2[i::3,j::3]
        Cmu[np.diag_indices_from(Cmu)] += sig_int**2 + self.dmz**2 
        
        C = inv(Cmu)
        
        L = self.distance_modulus_salt2(alpha, beta, Mb) - self.distance_modulus()
        self.residuals = L
        self.var = np.diag(Cmu)
        return P.dot(L,P.dot(C,L))





    def _compute_dispertion(self, sig_int):
        sig = optimize.fmin(self._disp_function, sig_int)[0]
        return sig
             
    def _disp_function(self,d):
        return abs((self.chi2_salt2(self.Params['alpha'], self.Params['beta'], self.Params['Mb'], d)/self.dof_salt2)-1.)
               
         
    def fit_sigma_int_salt2(self):
        
        sig_int = 0.001
        
        Find_param = minuit.Minuit(self.chi2_salt2, alpha=0.15,beta=2.9, Mb=-19.1, sig_int=sig_int, fix_sig_int= True)
        Find_param.migrad()
        self.Params = Find_param.values
        self.Params_Covariance = Find_param.covariance
        

        calls=0
        if abs((self.chi2_salt2(self.Params['alpha'], self.Params['beta'], self.Params['Mb'], sig_int)/(self.dof_salt2))-1.)>0.1:
            while abs((self.chi2_salt2(self.Params['alpha'], self.Params['beta'], self.Params['Mb'], sig_int)/(self.dof_salt2))-1.) > 0.001:
   
                if calls<100:
                    print 'je cherche de la dispersion pour la %i eme fois'%(calls+1)
                    sig_int = self._compute_dispertion(sig_int)                  
                    Find_param = minuit.Minuit(self.chi2_salt2, alpha=self.Params['alpha'], beta=self.Params['beta'], Mb=self.Params['Mb'], sig_int = sig_int, fix_sig_int= True)
                
                    Find_param.migrad()
                    self.Params = Find_param.values
                    print self.Params
                    calls+=1

                else:
                    print 'error : calls limit are exceeded'
                    break
        
        self.wrms_salt2, self.wrms_salt2_err = comp_rms(self.residuals, self.dof_salt2, err=True, variance=self.var)
        return Find_param, sig_int

    def _compute_dispertion_sugar(self, sig_int):
        sig = optimize.fmin(self._disp_function_sugar, sig_int)[0]
        print sig
        return sig
             
    def _disp_function_sugar(self,d):
        print d
        return abs((self.chi2_sugar(self.Params['cst'], self.Params['alpha1'], self.Params['alpha2'], self.Params['alpha3'], self.Params['beta'] ,d)/self.dof_sugar)-1.)        
        
        
    def fit_sigma_int_sugar(self):
        sig_int = 0.001
        
        if self.qi == False: 
            Find_param = minuit.Minuit(self.chi2_sugar, cst=0.01,alpha1=0., alpha2=0., alpha3=0., beta=0.,sig_int=sig_int, fix_sig_int= True, fix_alpha1=True, fix_alpha2=True, fix_alpha3=True, fix_beta=True)
            print sig_int
            Find_param.migrad()
            self.Params = Find_param.values
            self.Params_Covariance = Find_param.covariance
            
    
            calls=0
            if abs((self.chi2_sugar(self.Params['cst'], 0., 0., 0., 0.,sig_int)/(self.dof_sugar))-1.)>0.1:
                while abs((self.chi2_sugar(self.Params['cst'],0., 0., 0., 0.,sig_int)/(self.dof_sugar))-1.) > 0.001:
       
                    if calls<100:
                        print 'je cherche de la dispersion pour la %i eme fois'%(calls+1)
                        sig_int = self._compute_dispertion_sugar(sig_int)
                        print sig_int
                        Find_param = minuit.Minuit(self.chi2_sugar, cst=self.Params['cst'], alpha1=0., alpha2=0., alpha3=0., beta=0., sig_int = sig_int, fix_sig_int=True, fix_alpha1=True, fix_alpha2=True, fix_alpha3=True, fix_beta=True)             
                        Find_param.migrad()
                        self.Params = Find_param.values
                        print self.Params
                        calls+=1
                        
    
                    else:
                        print 'error : calls limit are exceeded'
                        break
        else:
            Find_param = minuit.Minuit(self.chi2_sugar, cst=0.01 ,alpha1=0.001, alpha2=0.001, alpha3=0.001, beta=0.001, sig_int=sig_int, fix_sig_int= True)
            print sig_int
            Find_param.migrad()
            self.Params = Find_param.values
            self.Params_Covariance = Find_param.covariance
            
    
            calls=0
            if abs((self.chi2_sugar(self.Params['cst'], self.Params['alpha1'], self.Params['alpha2'],self.Params['alpha3'], self.Params['beta'],sig_int)/(self.dof_sugar))-1.)>0.1:
                while abs((self.chi2_sugar(self.Params['cst'], self.Params['alpha1'], self.Params['alpha2'],self.Params['alpha3'], self.Params['beta'], sig_int)/(self.dof_sugar))-1.) > 0.001:
       
                    if calls<100:
                        print 'je cherche de la dispersion pour la %i eme fois'%(calls+1)
                        sig_int = self._compute_dispertion_sugar(sig_int)
                        print sig_int
                        Find_param = minuit.Minuit(self.chi2_sugar, cst=self.Params['cst'], alpha1=self.Params['alpha1'], alpha2=self.Params['alpha2'], alpha3=self.Params['alpha3'], beta=self.Params['beta'], sig_int = sig_int, fix_sig_int= True)             
                        Find_param.migrad()
                        self.Params = Find_param.values
                        print self.Params
                        calls+=1
                        
    
                    else:
                        print 'error : calls limit are exceeded'
                        break
            
        self.wrms_sugar, self.wrms_sugar_err = comp_rms(self.residuals, self.dof_sugar, err=True, variance=self.var)        
        return Find_param, sig_int
        

    def sig_fix_sugar(self, sig_int):
        if self.qi:
            Find_param = minuit.Minuit(self.chi2_sugar, cst=0.01 ,alpha1=0.001, alpha2=0.001, alpha3=0.001, beta=0.001, sig_int=sig_int, fix_sig_int= True)
        else:
            Find_param = minuit.Minuit(self.chi2_sugar, cst=0.01,alpha1=0., alpha2=0., alpha3=0., beta=0.,sig_int=sig_int, fix_sig_int= True, fix_alpha1=True, fix_alpha2=True, fix_alpha3=True, fix_beta=True)
        print sig_int
        Find_param.migrad()
        self.Params = Find_param.values
        self.Params_Covariance = Find_param.covariance       
        self.wrms_sugar, self.wrms_sugar_err = comp_rms(self.residuals, self.dof_sugar, err=True, variance=self.var)   
        return Find_param
    
    def plot_sig_wc_sugar(self):
        list_sig = np.linspace(0.09, 0.30, 30)
        res_wrms = []
        res_wrms_err = []
        res_chi2 = []
        for j in range(len(list_sig)):
            res = self.sig_fix_sugar(list_sig[j])
            res_wrms.append(self.wrms_sugar)
            res_wrms_err.append(self.wrms_sugar_err)
            res_chi2.append(res.fval/self.dof_sugar)
        
        p1 = plt.errorbar(list_sig, res_wrms, yerr=res_wrms_err, fmt='.')
        plt.xlabel( 'Sig int')
        plt.ylabel('wrms')
        ax2 = plt.gca().twinx()
        p2 = ax2.scatter(list_sig, res_chi2, color='red', marker='.')
        ax2.set_ylabel('chi2')
        if self.qi:
            pdffile = '../sugar_analysis_data/results/wrms_sig.pdf'
        else:
            pdffile = '../sugar_analysis_data/results/wrms_sig_Mgr.pdf'    
        plt.savefig(pdffile, bbox_inches='tight')  
        plt.show()

    def plot_HD(self,zcmb,dMgr ,dq1, dq2, dq3, dA):
        #plot setup
        rcParams['font.size'] = 16.
        font = {'family': 'normal', 'size': 16}
        rc('axes', linewidth=1.5)
        rc("text", usetex=True)
        rc('font', family='serif')
        rc('font', serif='Times')
        rc('legend', fontsize=25)
        rc('xtick.major', size=5, width=1.5)
        rc('ytick.major', size=5, width=1.5)
        rc('xtick.minor', size=3, width=1)
        rc('ytick.minor', size=3, width=1)
        
        mu = self.distance_modulus_sugar_corr(self.Params['cst'], self.Params['alpha1'], self.Params['alpha2'], self.Params['alpha3'], self.Params['beta'])
        dmu = self.distance_modulus_sugar_corr_uncertainty(self.Params['alpha1'], self.Params['alpha2'], self.Params['alpha3'], self.Params['beta'],dMgr ,dq1, dq2, dq3, dA)
        plt.errorbar(zcmb, mu, yerr=dmu, fmt='.')
        z = np.linspace(0.,0.12,200)
        self.zcmb = z
        zhl = self.zhl
        self.zhl = z
        plt.plot(z, self.distance_modulus(), color='r')
        plt.xlim(0.01,0.12)
        plt.ylim(32,40)
        plt.xlabel('z')
        plt.ylabel('$\mu = M_{gr} - cst - \sum_{i=1}^{3}  a_{i}q_{i} - b A_{v} $')
        self.zcmb = zcmb
        self.zhl = zhl
        pdffile = '../sugar_analysis_data/results/HD_sug.pdf'
        plt.savefig(pdffile, bbox_inches='tight')
         
