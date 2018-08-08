#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 13:02:35 2018

@author: florian
"""

import numpy as np
from operator import itemgetter
import sncosmo
from sugar_analysis import builtins_SNF as Build_SNF
from sncosmo.salt2utils import BicubicInterpolator
from scipy.interpolate import InterpolatedUnivariateSpline as Spline1d
from matplotlib import pyplot as plt
import copy
from astropy.table import Table
from sugar_analysis import math_toolbox as sam

CLIGHT = 2.99792458e18         # [A/s]
HPLANCK = 6.62606896e-27       # [erg s]
# wavelength limits for sugar model
wl_min_sug = 3254.01639
wl_max_sug = 8649.03871
wl_min_salt = 2000.000000
wl_max_salt = 9200.000000
t_min_sug = -12
t_max_sug = 48


Build_SNF.register_SNf_bands_width(width=10)
Build_SNF.mag_sys_SNF_width(width=10)
Build_SNF.register_SUGAR()
sugar_model = '/home/florian/sugar_model/'
sugar_analysis_data = '/home/florian/sugar_analysis_data/'
class sugar_spectrum():
    
    def __init__(self, modeldir=None,
                 m0file='sugar_template_0.dat',
                 m1file='sugar_template_1.dat',
                    m2file='sugar_template_2.dat',
                    m3file='sugar_template_3.dat',
                 m4file='sugar_template_4.dat', 
                 parameters_init=np.array([0., 0., 0., 0., 37.]),
                 name=None, version=None):
        

        self._SCALE_FACTOR = 1.
        self._param_names = ['q1', 'q2', 'q3', 'A', 'Mgr']
        
    
        infile = open(sugar_model+ 'SUGAR_model_v1.asci', 'r')
        lines = infile.readlines()
        infile.close()
        
        listik = []
        for line in lines:
            new_line = [float(i) for i in line.split()]
            listik.append(new_line)
        
        s = sorted(listik, key=itemgetter(0))
        
        names = ['0','1','2','3','4']
        for i in names:
            outfile = open(sugar_model + 'sugar_template_' + i + '.dat', 'w')
            for line in s:
                j = 2+int(i)
                outfile.write('%4.4f %8.8f %8.8f' %(line[0],line[1],line[j]))
                outfile.write('\n')
            outfile.close()        
            
        
        self.step = 15
        self.name = name
        self.version = version
        self._model_raw = {}
        self._model = {}
        self._parameters = parameters_init
        self._salt_params = [1., 0., 0.]
        names_or_objs = {'M0': m0file, 'M1': m1file, 'M2': m2file, 'M3': m3file, 'M4': m4file}

        # model components are interpolated to 2nd order
        for key in ['M0', 'M1', 'M2', 'M3', 'M4']:
            phase, wave, values = sncosmo.read_griddata_ascii(sugar_model + names_or_objs[key])
            values *= self._SCALE_FACTOR
            self._model_raw[key] = np.loadtxt(sugar_model+names_or_objs[key])
            self._model[key] = BicubicInterpolator(phase, wave, values)

            # The "native" phases and wavelengths of the model are those
            # of the first model component.
            if key == 'M0':
                self._phase = phase
                self._wave = wave
        self.dic_sts = None
    def model_spectrum_flux(self, phase, wave):
        '''
        '''
        m0 = self._model['M0'](phase, wave)
        m1 = self._model['M1'](phase, wave)
        m2 = self._model['M2'](phase, wave)
        m3 = self._model['M3'](phase, wave)
        m4 = self._model['M4'](phase, wave)
#        return 10. ** (-0.4 * (m0 + 48.59)) / (wave** 2 / 299792458. * 1.e-10)
        return 10. ** (-0.4 * (m0 + self._parameters[0] * m1 + self._parameters[1] * m2 + self._parameters[2] * m3 + self._parameters[3] * m4 + self._parameters[4]  + 48.59))/ (wave** 2 / 299792458. * 1.e-10)


#        source_salt2 = sncosmo.get_source('salt2')
#        flux_salt2 = source_salt2._flux(phase, wave)
#        return flux_salt2  

    def spectrum_generator(self, parameters=None, phase=None, mean_errors=None):
        """
        simulate sugar spectrum given a set of parameters
        parameters : array of sugar parameters (np.array([q1, q2, q3, A, Mgr]))
        mean_error : array of errors mean for different phase 
        """
        self.dic_spectrum = {}
        par_init = copy.deepcopy(self._parameters)
        if  not isinstance(phase, np.ndarray):
            phase = self._phase
        if  isinstance(parameters, np.ndarray):    
            self._parameters = parameters
        wave = np.linspace(3000,9000,10000)
        flux = self.model_spectrum_flux(phase, wave)
        flux_err = np.zeros_like(flux)
        self.dic_spectrum['parameters'] = parameters
        self.dic_spectrum['wave'] = wave
        self.dic_spectrum['phase'] = phase
        if not isinstance(mean_errors,np.ndarray):
            self.dic_spectrum['fluxerr'] = None
            self.dic_spectrum['flux'] = flux
        else:
            for i in range(len(flux)):
                flux_err[i]=np.ones(len(flux[i]))*mean_errors[i]
                for j in range(len(flux[i])):
                    flux[i,j] = np.random.normal(flux[i,j],mean_errors[i],1)[0]
            self.dic_spectrum['fluxerr'] = flux_err 
            self.dic_spectrum['flux'] = flux            
        self._parameters = par_init 
        
    def error_generator_lc_factory(self, band):
        """
        """
        if band == 'new_fI_10' or band =='fI_10':
            error = np.abs(np.random.normal(0.0101, 0.0425, 1)[0])
        if band == 'fU_10':
            error = np.abs(np.random.normal(0.0146, 0.0741, 1)[0])
        if band == 'fB_10':
            error = np.abs(np.random.normal(0.0141, 0.06898, 1)[0])
        if band == 'fV_10':
            error = np.abs(np.random.normal(0.0139, 0.0646, 1)[0])
        if band == 'fR_10':
            error = np.abs(np.random.normal(0.0088, 0.0395, 1)[0])
        return error
    
    def AstropyTable_flux(self, error=None):
        """
        """
    
                
        band_used = ['new_fI_10','fB_10','fV_10','fR_10','fU_10']    
        time = []
        flux = []
        fluxerr = []
        band = []
        zp = []
        zpsys = []
        vega = sncosmo.get_magsystem('vega_snf_10')
        for i in range(len(self._phase)):
            for b in band_used:
            
                filt2 = np.genfromtxt(sugar_analysis_data+'data/Instruments/Florian/'+b+'.dat')
                wlen = filt2[:,0]
                tran = filt2[:,1]
                self.splB = Spline1d(wlen, tran, k=1,ext = 1)    
                
                #computation of the integral
                dt = 100000
                xs = np.linspace(float(wl_min_sug), float(wl_max_sug), dt)
                self.xs = xs
                dxs = (float(wl_max_sug-wl_min_sug)/(dt-1))
                if error=='spectra':
                    spec_flux = np.random.normal(self.model_spectrum_flux(self._phase, xs),1e-34)
                else:
                    spec_flux = self.model_spectrum_flux(self._phase, xs)
#                plt.plot(xs,self.splB(xs))
#                plt.plot(xs,spec_flux[i,:])
#                print spec_flux[i,:] 
                plt.show()
                inte = np.sum(spec_flux[i,:]*(xs  / (CLIGHT*HPLANCK))*self.splB(xs)*dxs)
                    
                flux.append(inte)
                if error=='factory':
                    fluxerr.append(self.error_generator_lc_factory(b))
                elif error=='spectra':
                    fluxerr.append(spec_err)
                else:
                    fluxerr.append(0.0000000001)
                time.append(self._phase[i])
                band.append(b)
                zp.append(2.5*np.log10(vega.zpbandflux(b)))
                zpsys.append('vega_snf_10')
        data = Table([time, band, flux, fluxerr, zp, zpsys], names=('time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'), meta={'name': 'data'})
        return data
                
    def fit_lc_sugar(self, data ):
        """
        """
        source = sncosmo.get_source('sugar')
        dust = sncosmo.CCM89Dust()
        model = sncosmo.Model(source=source, effects=[dust],effect_names=['mw'], effect_frames=['obs'])    
        model.set(mwebv=0)
        model.set(z=0)
                   
        #initial iteration x1 fix
        print 'initialisation'            

        model.set(q1=0.0)
        model.set(q2=0.0)
        model.set(q3=0.0)     
        model.set(t0=0.0)
        res, fitted_model = sncosmo.fit_lc(data, model, ['A','Mgr'], modelcov = False)

            
        print 'first iteration'
        chi2 = res.chisq
        print chi2
        chi2p = chi2*2
        m=0
    
        while chi2 < chi2p and m < 10:

            print m
            if m > 0:
                resp = res
                fitted_modelp = fitted_model
            m += 1
            t_peak = fitted_model.parameters[1]
            #print t_peak,fitted_model.parameters[4]

            t1 = t_peak + t_min_sug*(1 + model.get('z'))
            t2 = t_peak + t_max_sug*(1 + model.get('z'))
                        
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



            print "WARNING: Sugar don't have model covariance for the moment"
            modelcov = False
            model.set(q1=res.parameters[2])
            model.set(q2=res.parameters[3])
            model.set(q3=res.parameters[4])
            model.set(A=res.parameters[5])
            model.set(Mgr=res.parameters[6])
            #print res.parameters[6]
            model.set(t0=res.parameters[1])                        
            res, fitted_model = sncosmo.fit_lc(data_new, model, ['q1', 'q2', 'q3', 'A', 'Mgr'], modelcov  = modelcov)
            

#                        sncosmo.plot_lc(data_new, model=fitted_model, errors=res.errors)
                    
#                        print res.chisq
    #                    print res

            chi2p = chi2
            chi2 = res.chisq
            print chi2p, chi2
        #final results
        res = resp
        fitted_model = fitted_modelp
        return res



    def get_param(self, sed, sed_err):
        
        W = np.eye(len(sed_err))*1./sed_err
        A = np.array([self._model_raw['M1'][:,2],
                      self._model_raw['M2'][:,2],
                      self._model_raw['M3'][:,2], 
                      self._model_raw['M4'][:,2], 
                      np.ones_like(self._model_raw['M0'][:,2])]).T
        left = np.linalg.inv(np.dot(A.T,np.dot(W, A)))
        right = np.dot(A.T, np.dot(W, sed - self._model_raw['M0'][:,2]))
        h = np.dot(left, right)
        return h
    
    def fit_spec_sugar(self):   
        spec_mag = []
        spec_mag_err = []
        for p in (self._phase):
            spec_mag_p = sam.flbda2ABmag(self._wave, self.model_spectrum_flux(p, self._wave)[0])  
            spec_mag_p_err = np.ones_like(spec_mag_p)
            spec_mag += list(spec_mag_p)
            spec_mag_err += list(spec_mag_p_err)
        spec_mag = np.array(spec_mag)
        spec_mag_err = np.array(spec_mag_err)    
        h = self.get_param(spec_mag, spec_mag_err)
        return h
    
    def mc_generator(self):
        """
        """
        data = self.AstropyTable_flux(error=True)
        
        for i in range(len(data['flux'])):
            if data['band'][i] == 'new_fI_10' or data['band'][i] =='fI_10':
                data['flux'][i] = np.random.normal(data['flux'][i], data['fluxerr'][i], 1)[0]
            if data['band'][i] == 'fU_10':
                data['flux'][i] = np.random.normal(data['flux'][i], data['fluxerr'][i], 1)[0]
            if data['band'][i] == 'fB_10':
                data['flux'][i] = np.random.normal(data['flux'][i],data['fluxerr'][i], 1)[0]
            if data['band'][i] == 'fV_10':
                data['flux'][i] = np.random.normal(data['flux'][i], data['fluxerr'][i], 1)[0]
            if data['band'][i] == 'fR_10':
                data['flux'][i] = np.random.normal(data['flux'][i], data['fluxerr'][i], 1)[0]

        return data
        



        
        
                