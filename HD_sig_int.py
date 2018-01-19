# -*- coding: utf-8 -*-
"""
Created on Fri Nov  3 12:20:25 2017

@author: mondon
"""

import numpy as np
import pylab as P
from scipy import integrate
from numpy.linalg import inv
import iminuit as minuit
from scipy import optimize
import copy



###########
#constants
###########

clight = 299792.458
H0 = 0.000070
###############
#
#def make_method(obj):
#    """Decorator to make the function a method of *obj*.
#    In the current context::
#      @make_method(Axes)
#      def toto(ax, ...):
#          ...
#    makes *toto* a method of `Axes`, so that one can directly use::
#      ax.toto()
#    COPYRIGHT: from Yannick Copin
#    """
#
#    def decorate(f):
#        setattr(obj, f.__name__, f)
#        return f
#
#    return decorate
#
#
#
def comp_rms(residuals, dof, err=True, variance=None):
    """
    Compute the RMS or WRMS of a given distribution.

    :param 1D-array residuals: the residuals of the fit.
    :param int dof: the number of degree of freedom of the fit.
    :param bool err: return the error on the RMS (WRMS) if set to True.
    :param 1D-aray variance: variance of each point. If given,
                             return the weighted RMS (WRMS).

#    :return: rms or rms, rms_err
#    """
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


#
#
#
###########################
##                        #
##   Hubble Fit           #
##                        #
########################### 
#PARAM_NAME = np.asarray(['alpha','beta',"gamma","epsilon","delta"])
#
#def get_hubblefit(x, cov_x, zhl, zcmb, zerr, parameters=None):
#    """
#    Parameters
#    ----------
#    x: type
#        infor
#    
#    cov_x ....
#    
#    
#    parameters: list of int to specify if you want to remove some parameters in the fit in PARAM_NAME. 
#                By default is None and this take all parameters. If you don't want to put a correction use 
#                parameters=[]
#        example: 
#            for usual salt x1 c correction if you want only color correction, use parameters = [2]. 
#            Warning if you have more than 5 parameters, add parameters in PARAM_NAME
#    """
#    n_corr =  np.shape(x)[1]-1
#    
#    class hubble_fit_case(Hubble_fit):
#        freeparameters = ["Mb"]+PARAM_NAME[:n_corr].tolist()
#        
#    if parameters is None:
#        parameters = np.arange(n_corr)
#        
#    h = hubble_fit_case(x, cov_x, zhl, zcmb, zerr)
#    # Do we have fixed (set to 0) parameters
#    for i,param in enumerate(h.freeparameters[1:]):
#        print 's'
#        if i not in parameters:
#            print 'l'
#            setattr(h,"%s_guess"%param,0)
#            setattr(h,"%s_fixed"%param,True)
#    return h
#    
#
#      
#class Hubble_fit(object):
#    
#    def __new__(cls,*arg,**kwargs):
#        """ Upgrade of the New function to enable the
#        the _minuit_ black magic
#        """
#        obj = super(Hubble_fit,cls).__new__(cls)
#        
#        exec "@make_method(Hubble_fit)\n"+\
#             "def _minuit_chi2_(self,%s): \n"%(", ".join(obj.freeparameters))+\
#             "    parameters = %s \n"%(", ".join(obj.freeparameters))+\
#             "    return self.get_chi2(parameters)\n"
#
#
#        return obj
#        
#    def __init__(self,X, cov_x, zhl, zcmb, zerr):
#        self.variable = np.asarray(X)
#        self.cov = cov_x
#        self.zcmb = zcmb
#        self.zhl = zhl
#        self.zerr = zerr
#        self.dmz = 5/np.log(10) * np.sqrt(self.zerr**2 + 0.001**2) / self.zcmb
#        self.dof = len(X)-len(self.freeparameters)  
#    
#    # ------------------ #
#    #   Cosmology        #
#    # ------------------ #
#
#    def int_cosmo(self, z, Omega_M=0.3):   
#        """
#        """
#        return 1./np.sqrt(Omega_M*(1+z)**3+(1.-Omega_M))
#        
#    def luminosity_distance(self):
#        """
#        """
#        if type(self.zcmb)==np.ndarray:
#            integr = np.zeros_like(self.zcmb)
#            for i in range(len(self.zcmb)):
#                integr[i]=integrate.quad(self.int_cosmo, 0, self.zcmb[i])[0]
#        else:
#            integr = integrate.quad(self.int_cosmo, 0, self.zcmb)[0]
#    
#        return (1+self.zhl)*(clight/H0)*integr
# 
#    def distance_modulus_th(self):
#        """
#        """
#        return 5.*np.log(self.luminosity_distance())/np.log(10.)-5.
#        
#    # ------------------ #
#    #   Fit              #
#    # ------------------ #
#    def setup_guesses(self,**kwargs):
#        """ Defines the guesses, boundaries and fixed values
#        that will be passed to the given model.
#        For each variable `v` of the model (see freeparameters)
#        the following array will be defined and set to param_input:
#           * v_guess,
#           * v_boundaries,
#           * v_fixed.
#        Three arrays (self.paramguess, self.parambounds,self.paramfixed)
#        will be accessible that will point to the defined array.
#        Parameter
#        ---------
#        **kwargs the v_guess, v_boundaries and, v_fixed for as many
#        `v` (from the freeparameter list).
#        All the non-given `v` values will be filled either by pre-existing
#        values in the model or with: 0 for _guess, False for _fixed, and
#        [None,None] for _boundaries
#        Return
#        ------
#        Void, defines param_input (and consquently paramguess, parambounds and paramfixed)
#        """
#        def _test_it_(k,info):
#            param = k.split(info)[0]
#            if param not in self.freeparameters:
#                raise ValueError("Unknown parameter %s"%param)
#
#        self.param_input = {}
#        # -- what you hard coded
#        for name in self.freeparameters:
#            for info in ["_guess","_fixed","_boundaries"]:
#                if hasattr(self, name+info):
#                    self.param_input[name+info] = eval("self.%s"%(name+info))
#                    
#        # -- Then, whatever you gave
#        for k,v in kwargs.items():
#            if "_guess" in k:
#                _test_it_(k,"_guess")
#            elif "_fixed" in k:
#                _test_it_(k,"_fixed")
#            elif "_boundaries" in k:
#                _test_it_(k,"_boundaries")
#            else:
#                raise ValueError("I am not able to parse %s ; not _guess, _fixed nor _boundaries"%k)
#            self.param_input[k] = v
#
#        # -- Finally if no values have been set, let's do it
#        for name in self.freeparameters:
#            if name+"_guess" not in self.param_input.keys():
#                self.param_input[name+"_guess"] = 0
#            if name+"_fixed" not in self.param_input.keys():
#                self.param_input[name+"_fixed"] = False
#            if name+"_boundaries" not in self.param_input.keys():
#                self.param_input[name+"_boundaries"] = [None,None]
#            
#    def fit(self, fit_intrinsic=True,**kwargs):
#        """
#        How to use kwargs 
#        For each variable `v` of the model (see freeparameters)
#        the following array will be defined and set to param_input:
#           * v_guess,
#           * v_boundaries,
#           * v_fixed.
#        Three arrays (self.paramguess, self.parambounds,self.paramfixed)
#        will be accessible that will point to the defined array.
#        Parameter
#        ---------
#        **kwargs the v_guess, v_boundaries and, v_fixed for as many
#        `v` (from the freeparameter list).
#        All the non-given `v` values will be filled either by pre-existing
#        values in the model or with: 0 for _guess, False for _fixed, and
#        [None,None] for _boundaries
#        """
#        self._loopcount = 0
#        self.sig_int = 0.
#        self.setup_guesses(**kwargs)
#        
#        first_iter = self._fit_minuit_()
#        # - Intrinsic disposerion Fit?
#        if fit_intrinsic:
#            while (np.abs(self.chi2_per_dof - 1) > 0.001 and self._loopcount < 10):
#                self.sig_int =  self.fit_intrinsic(np.sqrt(np.mean(self.var))*2. / self.chi2_per_dof)
#                iter = self._fit_minuit_()
#                self._loopcount += 1
#                
#        # - Final steps      
#        self._fit_readout_()
#        return 
#      
#    def get_chi2(self, params):
#        """
#        """
#        self.Cmu = np.zeros_like(self.cov[::len(params),::len(params)])
#        pcorr = np.concatenate([[1],params[1:]])
#        for i, coef1 in enumerate(pcorr):
#            for j, coef2 in enumerate(pcorr):
#                self.Cmu += (coef1 * coef2) * self.cov[i::len(params),j::len(params)]  
#                
#                
#        self.Cmu[np.diag_indices_from(self.Cmu)] += self.sig_int**2 + self.dmz**2 
#        self.C = inv(self.Cmu)
#        L = self.distance_modulus(params) - self.distance_modulus_th()
#        self.residuals = L
#        self.var = np.diag(self.Cmu)
#        return P.dot(L,P.dot(self.C,L))
#        
#    def distance_modulus(self, params):
#        """
#        (mb + alpha * v1 + beta * v2 .....) - Mb
#        """
#        return  np.sum(np.concatenate([[1],params[1:]]).T * self.variable, axis=1) - params[0]
#        
#        
#    def fit_intrinsic(self, intrinsic_guess=0.1):
#        """ Get the most optimal intrinsic dispersion given the current fitted standardization parameters. 
#        
#        The optimal intrinsic magnitude dispersion is the value that has to be added in quadrature to 
#        the magnitude errors such that the chi2/dof is 1.
#        Returns
#        -------
#        float (intrinsic magnitude dispersion)
#        """
#        def get_intrinsic_chi2dof(intrinsic):
#            self.sig_int = intrinsic
#            return np.abs(self.get_chi2(self.resultsfit) / self.dof -1)
#        
#        return optimize.fmin(get_intrinsic_chi2dof,
#                             intrinsic_guess, disp=0)[0]           
#    # ------------------ #
#    # Internal Fit Tools #
#    # ------------------ #
#    def _fit_minuit_(self,verbose=False, step=1):
#        """
#        """
#        self._setup_minuit_(step=step)
#        if verbose: print("STARTS MINUIT FIT")
#        self._migrad_output_ = self.minuit.migrad()
#        
#        if self._migrad_output_[0]["is_valid"] is False:
#            print("migrad is not valid")
#            self.fitOk = False
#        elif verbose:
#            self.fitOk = True
#            
#        self.resultsfit = np.asarray([self.minuit.values[k]
#                              for k in self.freeparameters])
#        self.chi2_per_dof = self.minuit.fval/self.dof
#                        
#            
#    def _setup_minuit_(self, step=1, print_level=0):
#        """
#        """
#        # == Minuit Keys == #
#        minuit_kwargs = {}
#        for param in self.freeparameters:
#            minuit_kwargs[param]           = self.param_input["%s_guess"%param]
#            minuit_kwargs["limit_"+param]  = self.param_input["%s_boundaries"%param]
#            minuit_kwargs["fix_"+param]    = self.param_input["%s_fixed"%param]
#
#        self.minuit = minuit.Minuit(self._minuit_chi2_,
#                             print_level=print_level,errordef=step,
#                             **minuit_kwargs)
#                        
#    def _fit_readout_(self):
#        """ Computes the main numbers """
#        
#        return comp_rms(self.residuals, self.dof, err=True, variance=self.var)   
#        
#        
#        
#        
#   """
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
        
class Hubble_fit_sugar(Hubble_fit):

    def __init__(self, X, cov_x, zhl, zcmb, zerr, qi=True):
        self.sugar_par = X
        self.cov = cov_x
        print cov_x
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
#            cov_mat = np.zeros([len(self.sugar_par), len(self.sugar_par)])
#            cv = self.cov
#    
#            for i in range (len(self.sugar_par)):
#                cov_mat[i,i] = cv[i,0,0]
#            self.cov_mat = cov_mat
            self.cov_mat = self.cov
            print self.cov
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
#            self.build_cov_mat()
          
            Cmu =  copy.deepcopy(self.cov)
            
            Cmu[np.diag_indices_from(Cmu)] += sig_int**2 + self.dmz**2 
            
            L = self.distance_modulus_grey(cst) - self.distance_modulus_th()

        else:
            self.build_cov_mat()
            print 'hello'
            Cmu = np.zeros_like(self.cov_mat[::5,::5])
            for i, coef1 in enumerate([1., -alpha1,-alpha2,-alpha3, -beta]):
                for j, coef2 in enumerate([1., -alpha1,-alpha2,-alpha3, -beta]):
                    Cmu += (coef1 * coef2) * self.cov_mat[i::5,j::5]
            Cmu[np.diag_indices_from(Cmu)] += sig_int**2 + self.dmz**2 
            L = self.distance_modulus_corr(cst, alpha1, alpha2, alpha3, beta) - self.distance_modulus_th()
            
#        self.Cmu = Cmu
        C = inv(Cmu)
        self.residuals = L
        self.var = np.diag(Cmu)            
        return P.dot(L.T,P.dot(C,L))        
        








    def _compute_dispertion(self, sig_int):
        sig = optimize.fmin(self._disp_function, sig_int)[0]
        print sig
        return sig
             
    def _disp_function(self,d):
        return abs((self.chi2(self.Params['cst'], self.Params['alpha1'], self.Params['alpha2'], self.Params['alpha3'], self.Params['beta'] ,d)/self.dof)-1.)        
        
        
    def fit_sigma_int(self):
        sig_int = 0.001
        
        if self.qi == False: 
            Find_param = minuit.Minuit(self.chi2, cst=0.01,alpha1=0., alpha2=0., alpha3=0., beta=0.,sig_int=sig_int, fix_sig_int= True, fix_alpha1=True, fix_alpha2=True, fix_alpha3=True, fix_beta=True)
#            print sig_int
            Find_param.migrad()
            self.Params = Find_param.values
#            self.Params_Covariance = Find_param.covariance
            print self.Params['cst']
    
            calls=0
            if abs((self.chi2(self.Params['cst'], 0., 0., 0., 0.,sig_int)/(self.dof))-1.)>0.1:
                while abs((self.chi2(self.Params['cst'],0., 0., 0., 0.,sig_int)/(self.dof))-1.) > 0.001:
       
                    if calls<100:
                        print 'je cherche de la dispersion pour la %i eme fois'%(calls+1)
                        sig_int = self._compute_dispertion(sig_int)
                        print sig_int
                        Find_param = minuit.Minuit(self.chi2, cst=self.Params['cst'], alpha1=0., alpha2=0., alpha3=0., beta=0., sig_int=sig_int, fix_sig_int=True, fix_alpha1=True, fix_alpha2=True, fix_alpha3=True, fix_beta=True)             
                        Find_param.migrad()
                        self.Params = Find_param.values
                        print Find_param.errors.get('cst')
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
#   """     